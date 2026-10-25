#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
bash scripts/build_linux.sh
bash scripts/run_from_zero_linux.sh "$@"
