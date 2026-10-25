#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ -n "${1:-}" ]];then
  logfile="$1"
elif [[ -f solver_live_windows.log && -f solver_live.log ]];then
  if [[ solver_live_windows.log -nt solver_live.log ]];then logfile=solver_live_windows.log;else logfile=solver_live.log;fi
elif [[ -f solver_live_windows.log ]];then logfile=solver_live_windows.log
else logfile=solver_live.log;fi
echo "Reading $logfile (run solver in another terminal)"
python3 visualization/cylinder_live_dashboard_windows.py --log "$logfile" --port 8765
