#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
[[ "${MSYSTEM:-}" == UCRT64 ]] || { echo 'Start MSYS2 UCRT64' >&2; exit 2; }
ck="${1:-}"
[[ -n "$ck" && -f "$ck/PARA" ]] || { echo 'Usage: bash scripts/run_from_checkpoint_ucrt64.sh examples/checkpoints/INC35_VERIFIED' >&2; exit 2; }
[[ -f solver_windows.exe ]] || bash scripts/build_windows_ucrt64.sh
python3 scripts/configure_run.py --mode restart
backup_run_data
cp -a "$ck"/. RESTART/
echo "Experimental continuation from checkpoint $ck; active self-contact pairs are NOT fully restorable."
run_with_log solver_windows.exe solver_live_windows.log
