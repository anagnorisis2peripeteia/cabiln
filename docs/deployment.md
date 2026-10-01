# Deployment and recovery

The public [CABILN app](https://cabiln.onrender.com/) automatically deploys from
`master`. Check a release by its Git SHA and the `/ready` response. The repository's
`deploy/render.yaml` is a separate staging template with automatic deployment
disabled; it does not describe the live service's deployment trigger.

## Release procedure

1. Push the candidate to a branch and check its **Tests** and **Fuzzing** workflows.
   All jobs must pass for the exact code being promoted.
2. Retain the candidate artifacts and the previous passing release outside CI's
   retention window.
3. Fast-forward `master` to the verified candidate. This triggers the live deploy.
4. Confirm that `/ready` reports the expected SHA. Exercise a drawing, notation
   conversion, Build preview/apply, highlighting, Undo and an exported molecule.
5. Check the workflows triggered by the `master` push as well.

```bash
release_ref=$(git rev-parse HEAD)
gh run list --repo anagnorisis2peripeteia/pyPept --commit "$release_ref"
curl -fsS https://cabiln.onrender.com/ready
```

The [CI workflows](../.github/workflows/) provide current release evidence.
[Earlier launch reviews](history.md) describe past checkpoints and their limits.
A passing functional check does not measure concurrent-user capacity: measure
editing latency, peak memory and cancellation under representative hosted load
before increasing concurrency.

## Build and retain an image

`Dockerfile` builds one wheel and installs the exact runtime dependencies from
`requirements-production.txt` with CPython 3.11.15. Build tools are pinned
separately. The Python base tag's OS layers can change, so retain the tested image.

```bash
docker build --build-arg CABILN_RELEASE="$(git rev-parse HEAD)" -t cabiln-candidate .
docker run --rm --memory=2g --cpus=1 -p 127.0.0.1:8000:8000 cabiln-candidate
```

Run `python tools/release_smoke.py` in another terminal. CI uses this production
profile, extracts its wheel for the Linux browser suite, and retains the image,
wheel, dependency hashes, library/rule binding and canonical convention for
30 days. Artifacts can appear before the workflow finishes; check the final
result before using them.

Render builds the live Git deployment separately. Verify that deployed revision
and its assets; the retained CI image remains a tested recovery artifact. For
an image-based service, publish the tested image and select its immutable digest.

## Hosting configuration

The staging Blueprint requests one 1-CPU/2-GiB instance, `/ready` health checks,
one chemistry worker, a 30-second deadline and a 1024-MiB worker address-space
limit. Applying it provisions a separate service. Use the actual service settings
when checking production capacity or storage.

Render documents deployment triggers in its
[Blueprint reference](https://render.com/docs/blueprint-spec#autodeploytrigger).
Use `/ready` for [health checks](https://render.com/docs/health-checks);
`/health` reports only that the application responds.

| Setting | Container default or supported limit |
| --- | --- |
| `CABILN_ENV` | `production`; requires process execution |
| `CABILN_EXECUTION` | `process` |
| `CABILN_WORKERS` | `1`; maximum 4 |
| `CABILN_JOB_TIMEOUT_SECONDS` | `30`; maximum 120 |
| `CABILN_WORKER_MEMORY_MB` | `1024`; Linux address-space limit |
| `CABILN_MAX_REQUEST_BYTES` | `2097152` |
| `CABILN_MAX_RESPONSE_BYTES` | `8388608` |
| `CABILN_ENABLE_REGISTRATION` | `0` |
| `CABILN_RELEASE` | Release ID; falls back to `RENDER_GIT_COMMIT` |

Use the `cabiln` entry point. It disables raw URL access logs and uses one Uvicorn
worker. Production suppresses native chemistry output. Logs contain request ID,
route, status, elapsed time and release; internal failures add exception type
and code locations without submitted structures. Set the host's proxy/log
retention separately. The application stores no server-side document database.

Busy calculations return 503 with `Retry-After`; deadlines terminate the worker
and return 504. Disconnects terminate active calculations. The browser makes up
to two cancellable retries for brief overload, then displays a Retry control.
Registration writes are never retried automatically. Workers prepare palette
data before reporting readiness; changing the library rebuilds that cache.

The render cache's 32-MiB budget counts retained payloads. Native allocations
and the parent process need additional memory. See [runtime](runtime-execution.md)
for limits, worker replacement, readiness and logging.

## Administrative ingestion and backups

Keep public registration disabled. For an administrative instance, copy the
bundled `monomers.sdf` and companion `monomers.csv` to a persistent data directory.
Set these before startup:

```text
CABILN_ENABLE_REGISTRATION=1
CABILN_MONOMER_LIBRARY=/data/library/monomers.sdf
CABILN_LIBRARY_BACKUP_DIR=/data/backups
CABILN_REGISTRATION_USER=admin
CABILN_REGISTRATION_TOKEN=<secret of at least 32 characters>
```

Supply the secret through the host's secret settings and use HTTPS. Both the
page and each write authenticate; opening the page does not authorize a later
unauthenticated write. Cross-origin writes are rejected. Without a token, the
local CLI only enables registration for loopback clients and refuses a public
bind. Production requires the external library, backup location and token.

Render's [persistent disks](https://render.com/docs/disks) preserve only their
mounted path; the default filesystem is ephemeral. A disk is tied to one
instance. Library and backups must be under durable mounts. Copy snapshots to
independent storage regularly; a backup on the same disk does not cover disk loss.

Every successful ingestion first takes a locked, versioned ZIP snapshot of the
previous SDF and aliases. A backup failure prevents the library write. Take a
manual snapshot with:

```bash
python -m pyPept.library_snapshots snapshot \
  --library /data/library/monomers.sdf --directory /data/backups
```

To restore, stop **all** readers and writers, select the matching release/rules,
then run:

```bash
python -m pyPept.library_snapshots restore --offline \
  --library /data/library/monomers.sdf --directory /data/backups \
  --snapshot /data/backups/library-<timestamp>-<id>.zip
```

Restore validates hashes, reaction binding and parsed definitions before writing,
takes a pre-restore snapshot, and rolls back if replacement fails. The SDF and
CSV are individually atomic replacements, so the service must be stopped while
the pair changes. `--offline` is an operator assertion; it does not discover
remote readers. Restart and confirm `/ready`, palette detection and a saved
project before reopening traffic.

For a Git-based deployment, redeploy the previous verified commit. An image-based
service can select its retained image digest. Use a compatible library snapshot
and preserve the current snapshot before restoring an older one. A
different canonical/RDKit convention intentionally requires source verification
when opening projects; do not silently relabel their chemistry metadata.
