FROM python:3.11.15-slim-bookworm AS wheel
WORKDIR /build
COPY requirements-build.txt ./
RUN python -m pip install --no-cache-dir -r requirements-build.txt
COPY pyproject.toml README.md LICENSE LICENSE-ORIGINAL requirements*.txt ./
COPY src ./src
RUN python -m build --wheel --no-isolation --outdir /dist

FROM python:3.11.15-slim-bookworm AS runtime
WORKDIR /opt/cabiln
ARG CABILN_RELEASE=""
ENV CABILN_RELEASE=${CABILN_RELEASE} \
    CABILN_ENV=production \
    CABILN_EXECUTION=process \
    CABILN_WORKERS=1 \
    CABILN_JOB_TIMEOUT_SECONDS=30 \
    CABILN_WORKER_MEMORY_MB=1024 \
    CABILN_ENABLE_REGISTRATION=0 \
    HOST=0.0.0.0 \
    PORT=8000 \
    OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1
RUN apt-get update \
    && apt-get install -y --no-install-recommends libxrender1 libxext6 \
    && rm -rf /var/lib/apt/lists/*
COPY requirements-production.txt ./
RUN python -m pip install --no-cache-dir --only-binary=:all: \
    --report dependency-install.json -r requirements-production.txt
COPY --from=wheel /dist ./dist
RUN python -m pip install --no-cache-dir --no-deps dist/*.whl \
    && python -m pip check
COPY tools/release_manifest.py ./release_manifest.py
RUN python release_manifest.py --wheel-dir dist --output release.json \
    && useradd --uid 10001 --create-home cabiln
USER cabiln
EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --start-period=45s \
    CMD python -c "import os,urllib.request; urllib.request.urlopen('http://127.0.0.1:'+os.environ['PORT']+'/ready',timeout=3).read()"
CMD ["cabiln"]
