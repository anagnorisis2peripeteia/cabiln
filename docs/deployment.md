# Deployment and recovery

The production candidate is a read-only public builder. Administrative ingestion
still feeds the same library, tiles, attachment detection and recognition.
No account/workspace service is introduced. A public multiuser private-monomer
service would require a separate ownership model.

## Release artifact

`Dockerfile` builds one wheel and installs it with the exact dependency closure
in `requirements-production.txt`, using CPython 3.11.15. The general library
requirements remain compatible with supported Python versions. Build tools are
pinned separately. The base tag pins Python, but its OS layers can change:
retain the built image by digest instead of assuming a rebuild is identical.

```bash
docker build --build-arg CABILN_RELEASE=<commit> -t cabiln-candidate .
docker run --rm --memory=2g --cpus=1 -p 127.0.0.1:8000:8000 cabiln-candidate
```

In a second terminal, run `python tools/release_smoke.py`. CI builds the image,
tests its running production profile, extracts that exact wheel for the Linux
browser suite, and retains the image, wheel, dependency-download hashes,
library/rule binding, canonical convention and image ID for 30 days. Artifacts
exist before all jobs finish, so an artifact is eligible for promotion only
when **all** jobs for its commit pass. Preserve the previous passing release
outside CI's retention window before promoting another release.

The frozen dependencies came from the tested local environment. A fresh Linux
binary-wheel installation and the container/browser jobs are mandatory gates;
see [current evidence](launch-validation.md). No Linux build is claimed merely
because versions are pinned.

## Hosting profile

`deploy/render.yaml` is a reviewable, manual-deploy staging Blueprint. It selects
one 1-CPU/2-GiB instance, `/ready`, one chemistry worker, a 30-second deadline,
and a 1024-MiB worker address-space limit. These are conservative starting
limits, not measured hosted capacity. Applying the Blueprint provisions paid
hosting. This implementation did not apply it or deploy the application.

Render's [Blueprint fields](https://render.com/docs/blueprint-spec) define the
plan, health route and disabled automatic deployment. For promotion, publish
the passing image and select its immutable registry digest in an image-backed
service. A Git-based Docker rebuild creates another candidate and needs staging
validation again. Render's [health checks](https://render.com/docs/health-checks)
should use `/ready`; `/health` deliberately remains cheap liveness.

| Setting | Production default or policy |
| --- | --- |
| `CABILN_ENV` | `production`; refuses unbounded local execution |
| `CABILN_EXECUTION` | `process` |
| `CABILN_WORKERS` | `1`; fixed pool, maximum 4 |
| `CABILN_JOB_TIMEOUT_SECONDS` | `30`; maximum 120 |
| `CABILN_WORKER_MEMORY_MB` | `1024`; Linux address-space bound, not RSS |
| `CABILN_MAX_REQUEST_BYTES` | `2097152` |
| `CABILN_MAX_RESPONSE_BYTES` | `8388608` |
| `CABILN_ENABLE_REGISTRATION` | `0` |
| `CABILN_RELEASE` | Commit/release ID; falls back to `RENDER_GIT_COMMIT` |

Use the `cabiln` entry point. It disables raw URL access logging and fixes the
parent to one Uvicorn worker. Native chemistry output is suppressed in production.
Logs contain request ID, route pattern, status, elapsed time and release. Internal
failures add exception type and code locations, without request bodies, molecule
text or exception messages. The app creates no server document database; browser
drafts and downloaded files have separate lifetimes. The host's own proxy/log
retention policy still needs to be set by the operator.

Busy work receives 503 with `Retry-After`; deadlines terminate the worker and
return 504. A disconnected request also terminates its calculation. Workers are
replaced after crashes and periodically recycled. The browser makes at most two
abort-aware retries for transient busy computation requests. Persistent overload
remains visible. It never automatically retries a registration write.

The 32-MiB render-cache budget counts retained Python payload size; it is not a
bound on total process memory. Linux enforces the worker address-space limit.
macOS tests exercise process lifecycle but do not prove that limit. Keep the
container memory cap, and measure total memory and normal editing latency under
representative concurrent work before increasing concurrency.

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
No storage was provisioned or inspected in this implementation.

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

For application rollback, select the retained previous image digest and its
compatible library snapshot. Preserve both current and previous snapshots. A
different canonical/RDKit convention intentionally requires source verification
when opening projects; do not silently relabel their chemistry metadata.
