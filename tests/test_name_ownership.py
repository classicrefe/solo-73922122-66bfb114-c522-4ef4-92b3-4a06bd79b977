"""Regression tests: lane names own exactly one directory, never more.

Covers the reported failures:
  * "." / ".." (and other non-component names) must be rejected before any
    upload is accepted, and must never touch other published lanes;
  * ".tmp.<x>" is a legitimate lane name and stays downloadable, includable
    in snapshots and safe from startup recovery;
  * genuinely interrupted staging is cleaned and a later restart/retry can
    continue;
  * a failing new task cannot desynchronise another lane's published row
    from its files.
"""
import gzip
import io
import json
import time

import pytest

BARCODES = [
    {"sample_id": "s1", "barcode_r1": "AAAAAAAA", "barcode_r2": "CCCCCCCC"},
]


def _post(client, lane, barcodes, r1_path, r2_path):
    with open(r1_path, "rb") as f1, open(r2_path, "rb") as f2:
        return client.post(
            "/jobs",
            data={"lane": lane, "barcodes": json.dumps(barcodes)},
            files={"r1": ("R1.fastq.gz", f1, "application/gzip"),
                   "r2": ("R2.fastq.gz", f2, "application/gzip")},
        )


def _wait(client, job_id, timeout=15):
    deadline = time.time() + timeout
    while time.time() < deadline:
        row = client.get(f"/jobs/{job_id}").json()
        if row["status"] in ("published", "failed"):
            return row
        time.sleep(0.1)
    raise AssertionError("timeout waiting for job")


def _pair(tmp):
    r1, r2 = tmp / "R1.fastq.gz", tmp / "R2.fastq.gz"
    with gzip.open(r1, "wb") as g:
        g.write(b"@x\nAAAAAAAATT\n+\nIIIIIIIIII\n")
    with gzip.open(r2, "wb") as g:
        g.write(b"@x\nCCCCCCCCGG\n+\nJJJJJJJJJJ\n")
    return r1, r2


def _publish(client, tmp, lane):
    r1, r2 = _pair(tmp)
    resp = _post(client, lane, BARCODES, r1, r2)
    assert resp.status_code == 201, resp.text
    done = _wait(client, resp.json()["job"]["job_id"])
    assert done["status"] == "published", done
    return done


def _read_gz(resp):
    with gzip.open(io.BytesIO(resp.content)) as g:
        return g.read()


def test_escaping_lane_names_rejected_before_touching_data(
        client, data_dir, tmp_path):
    """A "."/".." submission must not unpublish or delete other lanes."""
    _publish(client, tmp_path, "keepLane")
    conn_path = data_dir / "published" / "lanes" / "keepLane"
    files_before = {
        p.name: p.read_bytes() for p in conn_path.iterdir() if p.is_file()}
    downloaded = {
        "r1": _read_gz(client.get("/lanes/keepLane/samples/s1/r1")),
        "r2": _read_gz(client.get("/lanes/keepLane/samples/s1/r2")),
    }

    r1, r2 = _pair(tmp_path)
    for bad_lane in (".", "..", "", "a/b", "../keepLane", "/abs"):
        resp = _post(client, bad_lane, BARCODES, r1, r2)
        # Empty form fields are rejected by request validation (422); the
        # rest reach the lane-name check (400).  All are refused pre-accept.
        assert resp.status_code in (400, 422), (bad_lane, resp.text)

    # Nothing was published or deleted by the rejected submissions.
    assert client.get("/lanes/keepLane").status_code == 200
    assert {p.name: p.read_bytes() for p in conn_path.iterdir()
            if p.is_file()} == files_before
    assert _read_gz(client.get("/lanes/keepLane/samples/s1/r1")) == downloaded["r1"]
    assert _read_gz(client.get("/lanes/keepLane/samples/s1/r2")) == downloaded["r2"]

    jobs = client.get("/jobs").json()
    assert [j["lane"] for j in jobs] == ["keepLane"]
    # No upload directory was created for any rejected name.
    incoming = data_dir / "incoming"
    assert [p.name for p in incoming.iterdir()] == []

    snap = client.post("/snapshots").json()
    assert snap["lanes"] == ["keepLane"]


def test_dot_sample_id_rejected(client, tmp_path):
    r1, r2 = _pair(tmp_path)
    bad = [{"sample_id": "..", "barcode_r1": "AAAAAAAA",
            "barcode_r2": "CCCCCCCC"}]
    resp = _post(client, "sampleDot", bad, r1, r2)
    assert resp.status_code == 400


def test_tmp_prefixed_lane_is_fully_usable_and_survives_recovery(
        client, data_dir, tmp_path):
    """'.tmp.instrument-42' is a legal lane, not a staging marker."""
    lane = ".tmp.instrument-42"
    _publish(client, tmp_path, lane)

    assert client.get(f"/lanes/{lane}").status_code == 200
    assert client.get(f"/lanes/{lane}/samples/s1/r1").status_code == 200
    assert client.get(f"/lanes/{lane}/quarantine/r2").status_code == 200
    lane_dir = data_dir / "published" / "lanes" / lane
    assert (lane_dir / "s1.r1.fastq.gz").is_file()

    snap = client.post("/snapshots").json()
    assert lane in snap["lanes"]
    assert client.get(
        f"/snapshots/{snap['snapshot_id']}/samples/s1/r1").status_code == 200

    # A second instance's startup recovery must keep it: it has a publish row.
    from app import processing
    notes = processing.recover_on_startup()
    assert lane_dir.exists() and (lane_dir / "s1.r1.fastq.gz").is_file()
    assert not any(lane in n for n in notes)
    assert client.get(f"/lanes/{lane}").status_code == 200


def test_interrupted_staging_cleanup_is_restart_safe(data_dir, helpers):
    """Only genuine staging debris is removed; cleanup is idempotent."""
    import importlib
    from app import config as cfg, db as dbmod, processing

    importlib.reload(cfg)
    importlib.reload(dbmod)
    dbmod.init_db()

    # A published lane whose *name* starts with ".tmp." must be preserved.
    lane_dir = cfg.LANES_DIR / ".tmp.published"
    lane_dir.mkdir(parents=True)
    (lane_dir / "s1.r1.fastq.gz").write_bytes(b"published")
    dbmod.get_conn().execute(
        "INSERT INTO lane_publishes (lane, job_id, published_at, seq) "
        "VALUES ('.tmp.published','oldjob',?,1)", (dbmod.utcnow(),))

    # Debris from an aborted publish: staging directory, no job row at all.
    orphan_staging = cfg.STAGING_DIR / "orphanjob"
    orphan_staging.mkdir(parents=True)
    (orphan_staging / "partial.gz").write_bytes(b"junk")

    # Debris from a truly interrupted job: processing status, free lock.
    fixture_dir = data_dir.parent / "fastq"
    fixture_dir.mkdir(exist_ok=True)
    r1, r2 = _pair(fixture_dir)
    upload = cfg.INCOMING_DIR / "crashupload"
    upload.mkdir(parents=True)
    up1, up2 = upload / "R1.fastq.gz", upload / "R2.fastq.gz"
    up1.write_bytes(r1.read_bytes())
    up2.write_bytes(r2.read_bytes())
    dbmod.create_job(
        {"job_id": "crashjob", "lane": "crashLane", "input_hash": "x",
         "barcode_table": "[]", "r1_incoming": str(up1),
         "r2_incoming": str(up2), "created_at": dbmod.utcnow()},
        BARCODES,
    )
    dbmod.get_conn().execute(
        "UPDATE jobs SET status='processing', instance_id='ghost', "
        "heartbeat_at=? WHERE job_id='crashjob'",
        ("2000-01-01T00:00:00+00:00",))
    work = cfg.WORK_DIR / "crashjob"
    work.mkdir(parents=True)
    (work / "lock").write_text("ghost\n")
    (work / "output").mkdir()
    crash_staging = cfg.STAGING_DIR / "crashjob"
    crash_staging.mkdir(parents=True)
    (crash_staging / "half.gz").write_bytes(b"junk")

    notes = processing.recover_on_startup()
    assert not orphan_staging.exists()
    assert not crash_staging.exists()
    assert not work.exists()
    row = dbmod.get_conn().execute(
        "SELECT status FROM jobs WHERE job_id='crashjob'").fetchone()
    assert row["status"] == "failed"
    # A retry can now start from scratch: uploads are retained.
    assert up1.is_file() and up2.is_file()
    assert lane_dir.is_dir() and (lane_dir / "s1.r1.fastq.gz").is_file()
    assert row is not None

    assert processing.recover_on_startup() == []
    assert lane_dir.is_dir()


def test_new_task_failure_keeps_other_lane_published(client, data_dir, tmp_path):
    """A failing job must not detach another lane's row from its files."""
    _publish(client, tmp_path, "aliveLane")
    lane_dir = data_dir / "published" / "lanes" / "aliveLane"
    before = {p.name: p.read_bytes() for p in lane_dir.iterdir()}

    bad_r1, bad_r2 = _pair(tmp_path)
    bad_r1.write_bytes(bad_r1.read_bytes()[:-6])  # corrupt gzip
    resp = _post(client, "doomedLane", BARCODES, bad_r1, bad_r2)
    assert resp.status_code == 201
    failed = _wait(client, resp.json()["job"]["job_id"])
    assert failed["status"] == "failed"

    assert client.get("/lanes/aliveLane").status_code == 200
    assert client.get("/lanes/aliveLane/samples/s1/r1").status_code == 200
    assert {p.name: p.read_bytes() for p in lane_dir.iterdir()} == before
    # No directory was left behind for the failed lane.
    assert not (data_dir / "published" / "lanes" / "doomedLane").exists()
    assert not (data_dir / "published" / ".staging").joinpath(
        failed["job_id"]).exists()


def test_stale_attempt_cannot_publish(data_dir):
    """A requeued/failed job can never write publish metadata."""
    import importlib
    from app import config as cfg, db as dbmod

    importlib.reload(cfg)
    importlib.reload(dbmod)
    dbmod.init_db()
    dbmod.create_job(
        {"job_id": "stalejob", "lane": "staleLane", "input_hash": "x",
         "barcode_table": "[]",
         "r1_incoming": str(cfg.INCOMING_DIR / "x" / "R1.fastq.gz"),
         "r2_incoming": str(cfg.INCOMING_DIR / "x" / "R2.fastq.gz"),
         "created_at": dbmod.utcnow()},
        BARCODES,
    )
    dbmod.get_conn().execute(
        "UPDATE jobs SET status='failed', instance_id=NULL WHERE job_id='stalejob'")
    with pytest.raises(dbmod.StaleJobError):
        dbmod.publish_job("stalejob", "staleLane", [
            {"sample_id": "s1", "r1_file": "x", "r2_file": "y",
             "read_pairs": 0},
        ])
    assert dbmod.get_conn().execute(
        "SELECT 1 FROM lane_publishes WHERE lane='staleLane'").fetchone() is None
    assert dbmod.get_conn().execute(
        "SELECT status FROM jobs WHERE job_id='stalejob'").fetchone()[
            "status"] == "failed"


def test_recovery_ignores_unsafe_lane_row(data_dir, tmp_path):
    """A legacy/tampered row with lane '..' must not derive a filesystem path."""
    import importlib
    from app import config as cfg, db as dbmod, processing

    importlib.reload(cfg)
    importlib.reload(dbmod)
    dbmod.init_db()

    (cfg.LANES_DIR / "realLane").mkdir(parents=True)
    (cfg.LANES_DIR / "realLane" / "keep.txt").write_text("keep")
    dbmod.get_conn().execute(
        "INSERT INTO lane_publishes (lane, job_id, published_at, seq) "
        "VALUES ('realLane','oldjob',?,1)", (dbmod.utcnow(),))

    r1, r2 = _pair(tmp_path)
    upload = cfg.INCOMING_DIR / "unsafeupload"
    upload.mkdir(parents=True)
    up1, up2 = upload / "R1.fastq.gz", upload / "R2.fastq.gz"
    up1.write_bytes(r1.read_bytes())
    up2.write_bytes(r2.read_bytes())
    dbmod.create_job(
        {"job_id": "unsafejob", "lane": "..", "input_hash": "x",
         "barcode_table": "[]", "r1_incoming": str(up1),
         "r2_incoming": str(up2), "created_at": dbmod.utcnow()},
        BARCODES,
    )
    dbmod.get_conn().execute(
        "UPDATE jobs SET status='processing', instance_id='ghost' "
        "WHERE job_id='unsafejob'")
    work = cfg.WORK_DIR / "unsafejob"
    work.mkdir(parents=True)
    (work / "lock").write_text("ghost\n")

    notes = processing.recover_on_startup()
    assert any("unsafe" in n for n in notes), notes
    # The lanes root and its siblings were not removed.
    assert (cfg.LANES_DIR / "realLane" / "keep.txt").is_file()
    assert cfg.LANES_DIR.is_dir()
    assert cfg.PUBLISHED_DIR.is_dir()
