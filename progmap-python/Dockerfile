FROM python:3.11-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    TF_CPP_MIN_LOG_LEVEL=2

RUN apt-get update \
    && apt-get install --no-install-recommends -y libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/progmap
COPY requirements-linux-runtime.txt requirements-linux-cpu.txt ./
RUN python -m pip install --upgrade pip setuptools wheel \
    && python -m pip install -r requirements-linux-cpu.txt

COPY pyproject.toml README.md README_EN.md LICENSE MANIFEST.in ./
COPY src ./src
RUN python -m pip install --no-deps . \
    && useradd --create-home --uid 10001 progmap \
    && mkdir -p /data /results \
    && chown -R progmap:progmap /data /results

USER progmap
ENTRYPOINT ["progmap"]
CMD ["--help"]

