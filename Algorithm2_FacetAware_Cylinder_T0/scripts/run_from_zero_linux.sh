#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
[[ -f solver_linux ]] || bash scripts/build_linux.sh
python3 scripts/configure_run.py --mode zero "$@"
backup_run_data
[[ ! -f RESTART/PARA ]] || { echo 'Refusing: restart state unexpectedly present' >&2; exit 3; }
echo 'Algorithm2 fresh start: INCNUMBER=0 and TTIME=0 before increment 1; NO checkpoint loaded.'
run_with_log solver_linux solver_live.log
