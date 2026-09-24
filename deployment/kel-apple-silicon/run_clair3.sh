#!/bin/bash
set -euo pipefail
setup_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CONDA_PREFIX="$setup_dir/envs/clair3-arm64"
export PATH="$CONDA_PREFIX/bin:/opt/homebrew/bin:$PATH"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
exec "$CONDA_PREFIX/bin/run_clair3.sh" "$@"
