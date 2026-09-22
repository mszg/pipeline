#!/bin/bash
set -euo pipefail
setup_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$setup_dir/envs/kel-native/bin:/opt/homebrew/bin:$PATH"
export MPLCONFIGDIR="$setup_dir/cache/matplotlib"
export XDG_CACHE_HOME="$setup_dir/cache"
export XDG_CONFIG_HOME="$setup_dir/config"
export XDG_DATA_HOME="$setup_dir/data"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
mkdir -p "$MPLCONFIGDIR" "$XDG_CONFIG_HOME" "$XDG_DATA_HOME"
cd "$setup_dir/../workflow"
exec snakemake --cores 2 --scheduler greedy --rerun-incomplete "$@"
