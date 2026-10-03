# Production chemistry execution

Run the public service through `cabiln` with `CABILN_ENV=production` and
`CABILN_EXECUTION=process`. The CLI runs one ASGI server process, admits at most
16 HTTP connections/tasks, and disables raw access logs. Its chemistry pool
uses a separate persistent process for each configured job slot. Trusted local
use and ordinary endpoint unit tests retain the direct-call implementation.

The production profile starts conservatively:

| Setting | Value | Meaning |
| --- | --- | --- |
| `CABILN_WORKERS` | `1` | Concurrent chemistry jobs; there is no waiting job queue |
| `CABILN_JOB_TIMEOUT_SECONDS` | `30` | Wall time including transport and chemistry |
| `CABILN_WORKER_MEMORY_MB` | `1024` | Linux worker address-space limit, in MiB |
| `CABILN_MAX_REQUEST_BYTES` | `2097152` | Entire request body, including project wrapper |
| `CABILN_MAX_RESPONSE_BYTES` | `8388608` | Uncompressed worker response |
| `CABILN_CACHE_BYTES` | `33554432` | Retained render payloads per process |
| `CABILN_RELEASE` | release ID | Shared release marker; Render commit is the fallback |

Query strings are limited to 64 KiB. Existing Pydantic and chemistry validation
still applies, including source-length, MOL-size, image-dimension and recognition
search bounds. A request-body limit does not replace chemical validation.

Only the explicit read/edit routes in `web/execution.py:CHEMISTRY_ROUTES` enter
the pool. They include rendering, recognition, notation conversion, monomer
previews/lookups, bond editing, and project preparation/validation. Workers run
the existing ASGI application, so endpoint validation and response behavior stay
shared with local use. Registration writes stay in the authenticated parent;
they authenticate before any body is read, then use the same streaming body cap
and upload deadline before model parsing. They are never cancelled halfway
through a persistence operation.

An occupied pool returns `503` with `Retry-After: 1`. Oversized transport data
returns `413`. A job deadline returns `504` after its process has been killed
and reaped. Disconnects and service shutdown also terminate active work. A
replacement warms its library independently before becoming available. Crashes,
internal failures and recycling after 100 completed jobs replace the process;
failed replacements retry with a bounded backoff. No abandoned computation
continues in an executor thread. If termination itself fails, the slot stays
unavailable instead of being reused.

IPC uses length-bounded JSON on a private socket, never submitted pickle data.
Native worker stdout/stderr is discarded, independently of IPC. BLAS/OpenMP
thread counts are set to one before spawned children import their libraries.
Linux workers apply `RLIMIT_AS` after imports and before library warm-up. This
limits virtual address space, **not RSS or total container memory**. Positive
memory limits fail startup on unsupported platforms; local macOS execution can
use zero explicitly without claiming enforcement. Production requires a
positive limit. Keep one chemistry worker in an initial 2 GiB container, and
measure on the actual Linux deployment before increasing concurrency. The
container must also cover the web process, library data, transport buffers and
operating-system overhead.

The render LRU has both a 200-entry limit and a retained-byte limit. Accounting
walks the payload's Python containers, strings and scalars; cache bookkeeping and
native chemistry allocations remain outside that number. Oversized entries are
returned to the client without being cached. Worker recycling and the Linux
process limit bound allocations outside the render cache.

Each worker also retains up to 32 MiB / 64 entries of serialized assemblies.
This includes the product, unused-port graph and endpoint labels. Reads create
detached RDKit molecules; native working copies and Python cache bookkeeping
remain outside that byte count and within the process limit. Library definition
content and reaction/perception fingerprints bind every entry.

## Health, readiness and logs

`GET /health` returns exactly `{"status":"ok"}` from the event loop and makes
no chemistry call. `GET /ready` returns `200` only after startup has checked the
selected SDF and aliases, all reaction steps, cap rules, the shipped library
quality manifest and required browser assets, and at least one chemistry worker
is alive. The quality manifest need not match a custom external library; changed
definitions retain their unreviewed status. Busy workers remain ready;
replacing the only worker temporarily returns `503`. Probes reuse the startup
result rather than parsing molecules. A missing or malformed library can never
be ready. Fix broken startup configuration and restart; readiness does not
silently switch to the bundled library.

`/server_id` remains compatible but now returns the stable release ID. It does
not vary with PID or worker replacement. Frontend state must not reload merely
because a server instance changes.

Application responses get `X-Request-ID`. JSON logs contain the request ID,
method, known route pattern, status, release and duration. Unexpected failures
return a generic `500` with that ID and log only exception type and code-frame
locations. Submitted bodies, query strings, exception messages and locals are
excluded. Production also disables native RDKit logging in the parent so
administrative parsing cannot echo structures. Expected input failures retain
their useful HTTP validation messages. Keep raw URL access logging disabled
when supplying a different ASGI runner.

## Verification

On the development macOS host, one standalone run on 2026-09-29 measured startup
validation at 1.268 s and a warmed child at 2.758 s. Parent peak RSS was
159.6 MiB; the single child peaked at 179.7 MiB. These are local measurements
without a memory limit, not hosted-service capacity claims.

| Render | Elapsed | Response bytes | Retained payload bytes |
| --- | ---: | ---: | ---: |
| Glycine | 0.112 s | 6,306 | 11,635 |
| `K.[G(4,2).[ac(1,2)]]-A` | 0.045 s | 25,956 | 34,665 |
| Bundled semaglutide example | 1.399 s | 278,046 | 302,408 |
| Cached semaglutide | 0.011 s | 278,046 | 302,408 |

`tests/test_runtime.py` exercises real spawned work, timeout termination,
overload, cancellation, crash recovery, HTTP disconnect, active shutdown,
health during chemistry, byte bounds, readiness and safe logs. Its Linux-only
test checks `/proc` for the actual 1 GiB hard address-space limit, renders real
chemistry under that limit, and verifies that an insufficient limit cannot
start. The release browser job runs these checks against the installed wheel
and pinned Linux dependencies. The measurements above describe the dated macOS
run; use each commit's CI result for Linux verification and hosted measurements
for service capacity.
