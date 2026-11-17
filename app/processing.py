"""Core demultiplexing: streaming paired gzip FASTQ split with crash safety.

Memory profile: every record is handled one at a time using line-based IO;
no whole-file buffers are allocated, so a 512 MiB lane fits in a 128 MiB
container together with the Python runtime.
"""
from __future__ import annotations

import errno
import fcntl
import gzip
import os
import shutil
import time
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from . import config, db, naming

QUARANTINE = "__quarantine__"
_VALID_BASES = frozenset(b"ACGTN")


class FastqFormatError(Exception):
    """Raised when either mate is malformed; the whole lane then fails."""


class SizeExceededError(Exception):
    """Raised when decompressed input exceeds the 512 MiB budget."""


class InvalidNameError(Exception):
    """A job carries a lane/sample name that is not a safe path component."""


class AttemptSuperseded(Exception):
    """This attempt no longer owns its job; it must not publish."""
    def __init__(self, job_id: str):
        super().__init__(f"job {job_id} is no longer processing")
        self.job_id = job_id


@dataclass
class Record:
    name: bytes   # full header line including '@', without newline
    seq: bytes    # bases only
    plus: bytes   # header line without newline
    qual: bytes   # qualities only, raw Phred+33 bytes


def _read_record(gf: gzip.GzipFile, mate: str, index: int) -> Record | None:
    """Read and validate one FASTQ record. Returns None at clean EOF."""
    name = gf.readline()
    if name == b"":
        return None
    seq_line = gf.readline()
    plus_line = gf.readline()
    qual_line = gf.readline()
    if not name.startswith(b"@") or seq_line == b"" or not plus_line.startswith(b"+"):
        raise FastqFormatError(
            f"{mate}: record {index} is not a well-formed FASTQ record")
    if qual_line == b"" or seq_line.count(b"\n") == 0:
        raise FastqFormatError(f"{mate}: record {index} is truncated")

    seq = seq_line.rstrip(b"\r\n")
    qual = qual_line.rstrip(b"\r\n")
    name_token = name.rstrip(b"\r\n")
    if len(seq) == 0:
        raise FastqFormatError(f"{mate}: record {index} has an empty sequence")
    if any(b not in _VALID_BASES for b in seq):
        raise FastqFormatError(
            f"{mate}: record {index} contains characters outside A/C/G/T/N")
    if len(seq) != len(qual):
        raise FastqFormatError(
            f"{mate}: record {index} sequence/quality length mismatch")
    if any(q < 33 or q > 126 for q in qual):
        raise FastqFormatError(
            f"{mate}: record {index} has quality bytes outside Phred+33 range")
    if len(name_token) < 2:
        raise FastqFormatError(f"{mate}: record {index} has an empty read name")
    return Record(name_token, seq, plus_line.rstrip(b"\r\n"), qual)


def _paired_records(r1_path: Path, r2_path: Path,
                    size_budget: int) -> Iterator[tuple[Record, Record]]:
    """Yield validated R1/R2 pairs, enforcing name/order equality and budget."""
    used = 0
    with open(r1_path, "rb") as f1, open(r2_path, "rb") as f2:
        try:
            with gzip.GzipFile(fileobj=f1) as g1, gzip.GzipFile(fileobj=f2) as g2:
                index = 0
                while True:
                    try:
                        r1 = _read_record(g1, "R1", index)
                        r2 = _read_record(g2, "R2", index)
                    except (OSError, EOFError, zlib.error) as exc:
                        raise FastqFormatError(f"gzip stream corrupted: {exc}") from exc
                    if r1 is None and r2 is None:
                        return
                    if r1 is None or r2 is None:
                        raise FastqFormatError(
                            "R1 and R2 contain different numbers of records")
                    n1 = r1.name[1:].split()[0]
                    n2 = r2.name[1:].split()[0]
                    if n1 != n2:
                        raise FastqFormatError(
                            f"record {index}: R1/R2 read names differ "
                            f"({n1.decode('ascii', 'replace')} vs "
                            f"{n2.decode('ascii', 'replace')})")
                    used += (len(r1.name) + len(r1.seq) + len(r1.plus)
                             + len(r1.qual) + 8)
                    used += (len(r2.name) + len(r2.seq) + len(r2.plus)
                             + len(r2.qual) + 8)
                    if used > size_budget:
                        raise SizeExceededError(
                            f"decompressed total exceeds "
                            f"{config.MAX_TOTAL_UNCOMPRESSED} bytes")
                    yield r1, r2
                    index += 1
        except FastqFormatError:
            raise
        except (OSError, EOFError, zlib.error) as exc:
            raise FastqFormatError(f"gzip stream corrupted: {exc}") from exc


@dataclass
class SampleKey:
    sample_id: str
    combined: bytes  # r1 barcode + r2 barcode, 16 bases


def _best_sample(seq1: bytes, seq2: bytes,
                 keys: list[SampleKey]) -> SampleKey | None:
    """Return unique sample within <=1 substitution over both mates, else None."""
    if len(seq1) < config.BARCODE_LEN or len(seq2) < config.BARCODE_LEN:
        return None
    probe = seq1[:config.BARCODE_LEN] + seq2[:config.BARCODE_LEN]
    best: SampleKey | None = None
    best_dist = config.BARCODE_LEN * 2 + 1
    tied = False
    for key in keys:
        ref = key.combined
        dist = 0
        for a, b in zip(probe, ref):
            if a != b:
                dist += 1
                if dist > best_dist:
                    break
        if dist < best_dist:
            best_dist = dist
            best = key
            tied = False
        elif dist == best_dist:
            tied = True
    if best is None or best_dist > 1 or tied:
        return None
    return best


class _Writers:
    """One gzip writer per (bucket, mate); files opened eagerly."""

    def __init__(self, out_dir: Path, sample_ids: list[str]):
        self._files: dict[tuple[str, str], gzip.GzipFile] = {}
        out_dir.mkdir(parents=True, exist_ok=True)
        buckets = sample_ids + [QUARANTINE]
        for bucket in buckets:
            for mate in ("r1", "r2"):
                fh = open(out_dir / f"{bucket}.{mate}.fastq.gz", "wb")
                self._files[(bucket, mate)] = gzip.GzipFile(
                    filename="", mode="wb", fileobj=fh, compresslevel=6)
        self._raw = {k: v.fileobj for k, v in self._files.items()}

    def write(self, bucket: str, mate: str, rec: Record, trimmed: bool) -> None:
        gf = self._files[(bucket, mate)]
        gf.write(b"@" + rec.name[1:] + b"\n")
        if trimmed:
            gf.write(rec.seq[config.BARCODE_LEN:] + b"\n")
            gf.write(b"+\n")
            gf.write(rec.qual[config.BARCODE_LEN:] + b"\n")
        else:
            gf.write(rec.seq + b"\n")
            gf.write(b"+\n")
            gf.write(rec.qual + b"\n")

    def close(self) -> None:
        # GzipFile does not own an externally supplied fileobj: close the
        # gzip layer (flush footer + buffered data), then the raw file.
        for key, gf in self._files.items():
            gf.close()
            self._raw[key].close()


def _work_dir(job_id: str) -> Path:
    return config.WORK_DIR / job_id


def _lock_path(job_id: str) -> Path:
    return _work_dir(job_id) / "lock"


def lane_directory(lane: str) -> Path:
    """The one published directory that belongs to ``lane``.

    Only a single ordinary path component is allowed, so the result is
    always a direct child of the lanes root; ``.``/``..``/names with a
    separator raise instead of resolving to the root, its parent or another
    lane's directory.
    """
    if not naming.is_valid_lane_name(lane):
        raise InvalidNameError(f"invalid lane name: {lane!r}")
    candidate = config.LANES_DIR / lane
    if candidate.parent != config.LANES_DIR:
        raise InvalidNameError(f"lane {lane!r} escapes the lanes root")
    return candidate


def _staging_dir(job_id: str) -> Path:
    """Staging namespace is keyed by job id, never by lane name."""
    if not naming.is_safe_component(job_id):
        raise InvalidNameError(f"invalid job id: {job_id!r}")
    return config.STAGING_DIR / job_id


def _validated_sample_ids(sample_rows) -> list[str]:
    sample_ids = [r["sample_id"] for r in sample_rows]
    for sid in sample_ids:
        if not naming.is_valid_sample_id(sid):
            raise InvalidNameError(f"invalid sample_id: {sid!r}")
        if sid.startswith("__"):
            raise InvalidNameError(f"reserved sample_id: {sid!r}")
    return sample_ids


def _remove_if_unpublished(path: Path) -> None:
    """Remove ``path`` only when no publication references it.

    Callers pass a directory already proven to belong to the lane being
    processed; the publication check is the final guard against deleting a
    lane that another task published.
    """
    if not path.exists():
        return
    row = db.get_conn().execute(
        "SELECT 1 FROM lane_publishes WHERE lane=?", (path.name,)).fetchone()
    if row is None:
        shutil.rmtree(path, ignore_errors=True)


def acquire_job_lock(job_id: str):
    """Exclusive flock held for the whole processing/publish cycle."""
    wd = _work_dir(job_id)
    wd.mkdir(parents=True, exist_ok=True)
    fh = open(_lock_path(job_id), "w")
    try:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as exc:
        fh.close()
        if exc.errno in (errno.EACCES, errno.EAGAIN):
            raise RuntimeError(f"job {job_id} is locked by another instance") from exc
        raise
    fh.write(f"{os.getpid()}\n")
    fh.flush()
    return fh


def process_job(job_row, instance_id: str) -> None:
    """Full lifecycle: split -> stage -> atomic publish -> metadata."""
    job_id = job_row["job_id"]
    lane = job_row["lane"]
    final_publish = lane_directory(lane)
    staging = _staging_dir(job_id)
    sample_rows = db.get_samples(job_id)
    keys = [
        SampleKey(
            sample_id=r["sample_id"],
            combined=r["barcode_r1"].encode() + r["barcode_r2"].encode(),
        )
        for r in sample_rows
    ]
    sample_ids = _validated_sample_ids(sample_rows)
    r1_in = Path(job_row["r1_incoming"])
    r2_in = Path(job_row["r2_incoming"])

    lock = acquire_job_lock(job_id)
    work = _work_dir(job_id)
    staged = work / "output"
    try:
        # Remove leftovers from a previous crashed attempt of this same job.
        if staged.exists():
            shutil.rmtree(staged)
        # Staging leftovers from an interrupted run are cleaned at startup,
        # never by interpreting a user-supplied lane name.
        if staging.exists():
            shutil.rmtree(staging)
        # An unpublished directory for this lane can only be debris from
        # this or an earlier attempt; a published one is left untouched.
        _remove_if_unpublished(final_publish)

        writers = _Writers(staged, sample_ids)
        counts = {s: 0 for s in sample_ids}
        counts[QUARANTINE] = 0
        try:
            for pair_index, (r1, r2) in enumerate(
                    _paired_records(r1_in, r2_in,
                                    config.MAX_TOTAL_UNCOMPRESSED)):
                assigned = _best_sample(r1.seq, r2.seq, keys)
                if assigned is None:
                    bucket = QUARANTINE
                    writers.write(bucket, "r1", r1, trimmed=False)
                    writers.write(bucket, "r2", r2, trimmed=False)
                else:
                    bucket = assigned.sample_id
                    writers.write(bucket, "r1", r1, trimmed=True)
                    writers.write(bucket, "r2", r2, trimmed=True)
                counts[bucket] += 1
                if pair_index % config.HEARTBEAT_FLUSH_EVERY == 0:
                    db.heartbeat(job_id, instance_id)
                    # Stop promptly if this attempt lost ownership (the job
                    # was failed/recovered or re-claimed by a newer attempt),
                    # so it can never publish on top of newer work.
                    if not db.job_owned_by(job_id, instance_id):
                        raise AttemptSuperseded(job_id)
        finally:
            writers.close()

        # Every bucket (samples + quarantine), both mates, must exist and be
        # a readable gzip file, even when the bucket received no reads.
        for bucket in sample_ids + [QUARANTINE]:
            for mate in ("r1", "r2"):
                p = staged / f"{bucket}.{mate}.fastq.gz"
                if not p.exists() or p.stat().st_size == 0:
                    raise FastqFormatError(f"missing output file {p.name}")
                try:
                    with gzip.open(p, "rb") as verify:
                        while verify.read(65536):
                            pass
                except (OSError, EOFError, zlib.error) as exc:
                    raise FastqFormatError(
                        f"output {p.name} is not a valid gzip: {exc}") from exc

        # Stage inside the published filesystem, then atomically rename into
        # the final lane directory.  Staging and lane names live in separate
        # namespaces (published/.staging/<job> vs published/lanes/<lane>).
        staging.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(staged), str(staging))
        # Never replace a directory that is already published (a duplicate
        # submission reuses the task; only debris can reach this point).
        if final_publish.exists():
            if db.get_conn().execute(
                    "SELECT 1 FROM lane_publishes WHERE lane=?",
                    (lane,)).fetchone() is not None:
                raise RuntimeError(f"lane {lane} is already published")
            shutil.rmtree(final_publish)
        # The job must still be ours at the last moment before the files
        # become visible; publish_job re-checks inside its transaction.
        if not db.job_owned_by(job_id, instance_id):
            raise AttemptSuperseded(job_id)
        os.rename(staging, final_publish)

        count_rows = [
            {
                "sample_id": bucket,
                "r1_file": str(final_publish / f"{bucket}.r1.fastq.gz"),
                "r2_file": str(final_publish / f"{bucket}.r2.fastq.gz"),
                "read_pairs": n,
            }
            for bucket, n in counts.items()
        ]
        try:
            db.publish_job(job_id, lane, count_rows, instance_id)
        except db.StaleJobError:
            # A stale attempt lost its job row: the files it just moved are
            # not referenced by any publication.  Remove only those; other
            # lanes' published directories are never in scope.
            _remove_if_unpublished(final_publish)
            raise

        # Publish is durable; only now remove inputs and scratch space.
        for p in (r1_in, r2_in):
            try:
                p.unlink()
            except FileNotFoundError:
                pass
        # The upload directory is private to this job (plus its COMPLETE
        # marker), so it can go once the lane is visible.
        shutil.rmtree(r1_in.parent, ignore_errors=True)
        shutil.rmtree(work, ignore_errors=True)
    except AttemptSuperseded:
        # Leave the stage for the next attempt/recovery; never publish.
        raise
    finally:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        finally:
            lock.close()


# ---------------- startup recovery ----------------

def _lock_is_free(path: Path) -> bool:
    """True iff nobody currently holds flock on path."""
    if not path.exists():
        return True
    fh = open(path, "r")
    try:
        try:
            fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            if exc.errno in (errno.EACCES, errno.EAGAIN):
                return False
            raise
        fcntl.flock(fh.fileno(), fcntl.LOCK_UN)
        return True
    finally:
        fh.close()


def recover_on_startup() -> list[str]:
    """Clean debris from truly interrupted runs. Never touches live or
    published data: flock ownership protects live instances, DB rows protect
    published lanes."""
    notes: list[str] = []
    conn = db.get_conn()
    jobs_by_id = {r["job_id"]: r for r in conn.execute(
        "SELECT * FROM jobs").fetchall()}

    # 1. Jobs stuck in 'processing' whose owner lock is free (process dead).
    #    The flock is the authority: a live instance always holds it, so a
    #    free lock means the owner truly disappeared even if the heartbeat
    #    timestamp happens to look fresh.
    for job in conn.execute(
            "SELECT * FROM jobs WHERE status='processing'").fetchall():
        lock = _lock_path(job["job_id"])
        if not _lock_is_free(lock):
            continue  # another instance is actively working this job
        notes.append(f"recovering interrupted job {job['job_id']}")
        conn.execute(
            "UPDATE jobs SET status='failed', instance_id=NULL,"
            "error='interrupted (owner process disappeared)', updated_at=? "
            "WHERE job_id=?", (db.utcnow(), job["job_id"]))
        _remove_job_scratch(job, notes)

    # 2. Orphan work directories.  A directory whose job row is still
    #    'processing' belongs to a live (or not yet recovered) owner and is
    #    only removed via step 1.  A free lock plus a job row not processing
    #    means the attempt is over; a lock-free directory with no row at all
    #    is debris from an aborted publish.
    for wd in config.WORK_DIR.iterdir():
        if not wd.is_dir():
            continue
        job = jobs_by_id.get(wd.name)
        if job is not None and job["status"] == "processing":
            continue  # live owner or handled by step 1
        lock = wd / "lock"
        if lock.exists() and not _lock_is_free(lock):
            continue  # owned by a live instance
        if job is None and time.time() - wd.stat().st_mtime < 5:
            continue  # just created, likely a racing sibling
        # Remove the job's staging area together with its scratch space, so
        # a truly interrupted publish is reclaimed in one restart pass.
        shutil.rmtree(config.STAGING_DIR / wd.name, ignore_errors=True)
        notes.append(f"removing orphan work directory {wd.name}")
        shutil.rmtree(wd, ignore_errors=True)

    # 3. Unpublished lane directories left by a crashed publish.
    #    A directory is removable only if no publication row names it and no
    #    live instance is currently publishing that lane.  The name is never
    #    interpreted (".tmp.<x>" is a perfectly legal lane name), so files of
    #    a successful ".tmp.instrument-42" publication are protected by the
    #    DB row exactly like any other lane.
    lanes_root = config.LANES_DIR
    if lanes_root.exists():
        published = {r["lane"] for r in conn.execute(
            "SELECT lane FROM lane_publishes").fetchall()}
        for child in lanes_root.iterdir():
            if not child.is_dir() or child.name in published:
                continue
            guarded = False
            for job in jobs_by_id.values():
                if job["status"] != "processing" \
                        or _lock_is_free(_lock_path(job["job_id"])):
                    continue
                # A live owner publishing this lane, or (for a deployment
                # upgraded mid-publish) still staging under the legacy
                # ".tmp.<job>" name: protect it either way.
                if job["lane"] == child.name \
                        or child.name == f".tmp.{job['job_id']}":
                    guarded = True
                    break
            if not guarded:
                notes.append(
                    f"removing unpublished lane directory {child.name}")
                shutil.rmtree(child, ignore_errors=True)

    # 3b. Publish staging areas (separate namespace, keyed by job id).
    #     Only staging whose job is not being actively processed is cleaned;
    #     a truly interrupted staging area is removed and will be re-created
    #     by the next attempt (restart-safe, idempotent).
    if config.STAGING_DIR.exists():
        for child in config.STAGING_DIR.iterdir():
            if not child.is_dir():
                continue
            job = jobs_by_id.get(child.name)
            if job is not None and job["status"] == "processing" \
                    and not _lock_is_free(_lock_path(child.name)):
                continue  # live publish in progress
            notes.append(f"removing unfinished publish staging {child.name}")
            shutil.rmtree(child, ignore_errors=True)

    # 4. Incoming debris of jobs that no longer need their uploads.
    #    A directory without a COMPLETE marker is an upload still streaming
    #    (or a freshly aborted one younger than the grace period) and is
    #    never touched, so another instance's in-flight upload is safe.
    keep_incoming = {
        Path(r["r1_incoming"]).parent.name
        for r in conn.execute(
            "SELECT r1_incoming FROM jobs WHERE status IN "
            "('queued','processing','failed')").fetchall()
    }
    grace_seconds = 3600
    if config.INCOMING_DIR.exists():
        now = time.time()
        for d in config.INCOMING_DIR.iterdir():
            if not d.is_dir() or d.name in keep_incoming:
                continue
            if not (d / "COMPLETE").exists():
                # Aborted mid-upload: only reclaim well after any plausible
                # transfer could still be running.
                if now - d.stat().st_mtime < grace_seconds:
                    continue
            notes.append(f"removing stale incoming upload {d.name}")
            shutil.rmtree(d, ignore_errors=True)
    return notes


def _remove_job_scratch(job_row, notes: list[str]) -> None:
    job_id = job_row["job_id"]
    lane = job_row["lane"]
    wd = _work_dir(job_id)
    if wd.exists():
        shutil.rmtree(wd, ignore_errors=True)
    tmp = _staging_dir(job_id)
    if tmp.exists():
        shutil.rmtree(tmp, ignore_errors=True)
    try:
        final = lane_directory(lane)
    except InvalidNameError:
        # The row predates validation or was tampered with: never derive a
        # filesystem path from it.
        notes.append(f"skipped unsafe lane name for job {job_id!r}")
        return
    if final.exists() and db.get_conn().execute(
            "SELECT 1 FROM lane_publishes WHERE lane=?",
            (lane,)).fetchone() is None:
        shutil.rmtree(final, ignore_errors=True)
    notes.append(f"cleared scratch for job {job_id}")
