#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
[[ "${MSYSTEM:-}" == UCRT64 ]] || { echo 'Start MSYS2 UCRT64' >&2; exit 2; }
[[ -f solver_windows.exe ]] || bash scripts/build_windows_ucrt64.sh
python3 scripts/configure_run.py --mode zero "$@"
backup_run_data
[[ ! -f RESTART/PARA ]] || { echo 'Refusing: restart state unexpectedly present' >&2; exit 3; }
echo 'Algorithm1 fresh start: INCNUMBER=0 and TTIME=0 before increment 1; NO checkpoint loaded.'
run_with_log solver_windows.exe solver_live_windows.log
