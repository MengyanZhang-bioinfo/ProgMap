#!/usr/bin/env bash
#SBATCH --job-name=progmap
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=progmap-%j.log

set -Eeuo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: sbatch $0 /absolute/path/to/PANCANCER /absolute/path/to/results" >&2
  exit 2
fi

project_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "$project_dir/.venv/bin/activate"

exec progmap \
  --data-root "$1" \
  --output "$2" \
  --cancers all \
  --test ttest \
  --top-n all \
  --device auto \
  --threads "${SLURM_CPUS_PER_TASK:-8}"

