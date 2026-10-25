#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 visualization/export_fem_outputs.py --root . --output PARAVIEW_EXPORT "$@"
