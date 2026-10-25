#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
ck="${1:-}"
[[ -n "$ck" && -f "$ck/PARA" ]] || { echo 'Usage: bash scripts/run_from_checkpoint_linux.sh examples/checkpoints/INC35_VERIFIED' >&2; exit 2; }
[[ -f solver_linux ]] || bash scripts/build_linux.sh
python3 scripts/configure_run.py --mode restart
backup_run_data
cp -a "$ck"/. RESTART/
echo "Experimental continuation from checkpoint $ck; active self-contact pairs are NOT fully restorable."
run_with_log solver_linux solver_live.log
