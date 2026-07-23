#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 /absolute/path/to/PANCANCER /absolute/path/to/results" >&2
  exit 2
fi

data_root="$1"
output_root="$2"
threads="${PROGMAP_THREADS:-8}"
cancers="${PROGMAP_CANCERS:-all}"
test_method="${PROGMAP_TEST:-ttest}"
top_n="${PROGMAP_TOP_N:-all}"
device="${PROGMAP_DEVICE:-auto}"

exec progmap \
  --data-root "$data_root" \
  --output "$output_root" \
  --cancers "$cancers" \
  --test "$test_method" \
  --top-n "$top_n" \
  --device "$device" \
  --threads "$threads"

