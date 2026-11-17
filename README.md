# Sequencing Lane Demultiplexing API

FastAPI + SQLite service that splits mixed paired-end gzip FASTQ uploads into
per-sample files. Designed to run as **two container instances on one shared
volume** with a **128 MiB memory limit**.

## What it does

* Accepts a lane number, a dual-index barcode table (1–8 samples, two
  8-base barcodes each) and paired `R1`/`R2` gzip FASTQ files.
  Decompressed input across both mates is capped at **512 MiB**.
* Validates both mates: valid gzip, well-formed FASTQ, equal record counts,
  identical read names in the same order, matching sequence/quality lengths.
  Any corruption fails the **whole lane**.
* Assigns a pair to a sample when the combined R1+R2 barcode matches with
  **at most one substitution** and the candidate sample is **unique**.
  Everything else (no hit, >1 substitution, or a tie) goes to quarantine.
* Assigned reads have both 8-base barcodes and their corresponding quality
  characters trimmed; quarantined reads are kept intact.
* A lane becomes visible only when **all** per-sample paired files, the
  quarantine paired files and the count table exist (atomic publish:
  files staged, renamed, then a single DB transaction).
* Lane names are **one ordinary path component** (`[A-Za-z0-9_.-]{1,64}`,
  not `.`/`..`, no separators). Such names are rejected with `400` before
  any upload is spooled, so they can never resolve to the lanes root, its
  parent, or another lane's directory. Names that merely start with `.tmp.`
  (for example `.tmp.instrument-42`) are ordinary, fully supported lanes.
* Same lane + identical input (hash of barcode table and both files) reuses
  the existing task; same lane with different input is rejected with `409`.
* Failed tasks can be retried (`POST /jobs/{id}/retry`); uploads are kept
  until a lane is published.
* Snapshots freeze the published lane list at creation time. Sample/
  quarantine downloads through a snapshot only ever concatenate files from
  that fixed list — lanes published later never mix in.

## Running (Docker Compose, two instances on a shared volume)

```bash
docker compose up --build
# instance A on http://localhost:8001, instance B on http://localhost:8002
# Swagger demo UI: http://localhost:8001/docs  and  http://localhost:8002/docs
```

Both containers mount the named volume `demux-shared` at `/data`, are capped
at 128 MiB (`mem_limit`), and automatically recover on restart.

### Quick demo (through /docs)

```bash
python3 scripts/make_demo_fastq.py /tmp/demo          # writes R1/R2.fastq.gz
python3 -m pip install --user -r requirements.txt     # host-side, for the helper
```

1. Open <http://localhost:8001/docs>, expand `POST /jobs`, click *Try it out*.
2. Set `lane` to `demoLane` and `barcodes` to the JSON below.
3. Upload `/tmp/demo/R1.fastq.gz` for `r1` and `/tmp/demo/R2.fastq.gz` for `r2`,
   then *Execute*. The response carries `job.job_id` with `status: queued`.

```json
[{"sample_id":"sampleA","barcode_r1":"AAAAAAAA","barcode_r2":"CCCCCCCC"},
 {"sample_id":"sampleB","barcode_r1":"GGGGGGGG","barcode_r2":"TTTTTTTT"}]
```

Then, in the same UI (`GET /jobs/{job_id}` shows `published`):

```bash
curl http://localhost:8001/lanes/demoLane                 # counts
curl -OJ http://localhost:8001/lanes/demoLane/samples/sampleA/r1
curl -OJ http://localhost:8001/lanes/demoLane/quarantine/r1
curl -X POST http://localhost:8001/snapshots              # freeze lane list
curl -OJ http://localhost:8001/snapshots/<id>/samples/sampleA/r1
```

The same lane is served identically by instance B on port 8002. A lane named
`.tmp.instrument-42` works exactly like `demoLane` (it is a normal name), and
`"."`/`".."` are refused with `400`.

## API surface

| Method | Path | Purpose |
|---|---|---|
| POST | `/jobs` | Submit a lane (multipart: `lane`, `barcodes` JSON, `r1`, `r2`); invalid/escaping lane names → `400` |
| GET | `/jobs`, `/jobs/{id}` | Task status (`queued`/`processing`/`published`/`failed`) |
| POST | `/jobs/{id}/retry` | Requeue a failed task |
| GET | `/lanes`, `/lanes/{lane}` | Published lanes and per-bucket counts |
| GET | `/lanes/{lane}/samples/{sample_id}/{r1\|r2}` | Download a sample mate |
| GET | `/lanes/{lane}/quarantine/{r1\|r2}` | Download a quarantine mate |
| POST | `/snapshots` | Freeze the current published lane list |
| GET | `/snapshots`, `/snapshots/{id}` | List/inspect snapshots |
| GET | `/snapshots/{id}/samples/{sample_id}/{mate}` | Snapshot-scoped concatenated download |
| GET | `/snapshots/{id}/quarantine/{mate}` | Snapshot-scoped quarantine download |
| GET | `/health` | Liveness, instance id, current job, queue depth |

## Concurrency and crash safety

* **SQLite WAL** with `BEGIN IMMEDIATE` job claiming: exactly one instance
  ever picks up a given job; the other keeps serving requests.
* Each instance processes one lane at a time in a background thread, so a
  slow lane never blocks uploads, queries or downloads on the event loop.
* A job runs under an exclusive `flock` on its work directory. Startup
  recovery removes scratch/staging debris only when no lock is held, and
  never deletes directories that have a published-lane DB row — so another
  instance's in-flight work and all published data are preserved.
* Publish order: write/verify outputs in `work/<job>/output`, move to the
  dedicated staging namespace `published/.staging/<job>` (keyed by job id,
  never by lane name), atomic `rename` to `published/lanes/<lane>`, then
  insert metadata in one DB transaction. Inputs are deleted only after the
  commit, so a crash at any point leaves no half-visible lane.
* Filesystem ownership follows the database, not the name: recovery deletes
  a lane directory only when no `lane_publishes` row names it and no live
  instance holds that lane's job lock. A lane legitimately named
  `.tmp.instrument-42` is therefore protected exactly like any other lane,
  while staging debris under `published/.staging/` is reclaimed.
* A stale processing attempt (its job was requeued/failed elsewhere) checks
  job ownership before and inside the publish transaction and can never
  write `lane_publishes`/`lane_counts`. A failing new lane never changes
  another lane's published row or files.

## Memory

All I/O is streamed record-by-record (1 MiB network chunks, line-based gzip
reads, one open gzip writer per bucket). No file is ever buffered whole.
The automated streaming test processes hundreds of MiB of decompressed
data with RSS growth below 64 MiB.

## Layout

```
app/            FastAPI app, SQLite layer, streaming splitter, worker, naming
scripts/        demo FASTQ generator, two-instance integration check
tests/          pytest suite (happy path, corruption, retry, recovery,
                barcode rules, size budget, streaming memory, name ownership)
Dockerfile      python:3.11-slim, non-root user, /data volume, healthcheck
docker-compose.yml  two 128 MiB instances sharing one named volume
```

## Tests

```bash
python3 -m pip install -r requirements.txt pytest httpx
python3 -m pytest tests/ -q
python3 scripts/integration_two_instances.py   # real two-process check
```

Regression coverage for lane-name ownership lives in
`tests/test_name_ownership.py`: `"."`/`".."`/separator names rejected before
acceptance without touching published lanes or files; `.tmp.instrument-42`
published, downloadable, snapshotted and preserved by another instance's
recovery; interrupted staging reclaimed idempotently across restarts; a
failing new lane leaving another lane's `published` row and files intact;
and a superseded attempt unable to publish.
