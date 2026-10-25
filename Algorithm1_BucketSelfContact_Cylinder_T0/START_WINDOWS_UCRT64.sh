#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
[[ "${MSYSTEM:-}" == "UCRT64" ]] || { echo 'Open MSYS2 UCRT64' >&2; exit 2; }
bash scripts/build_windows_ucrt64.sh
bash scripts/run_from_zero_ucrt64.sh "$@"
