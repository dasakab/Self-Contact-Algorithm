#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

backup_run_data() {
  local stamp
  stamp="$(date -u +%Y%m%dT%H%M%SZ)_$$_$RANDOM"
  mkdir -p "run_backups/$stamp"
  for part in RESTART OUT solver_live.log solver_live_windows.log monitor.html contact_ring_events.csv PARAVIEW_EXPORT monitor_events.log ring_events.log monitor_events_windows.log ring_events_windows.log SECOND_RING_CANDIDATE.txt; do
    if [[ -e "$part" ]]; then
      mv "$part" "run_backups/$stamp/"
    fi
  done
  mkdir -p RESTART OUT PARAVIEW_EXPORT
  touch RESTART/.gitkeep OUT/.gitkeep
  echo "Previous outputs, if any, backed up under run_backups/$stamp/"
}

run_with_log() {
  local executable="$1"; local log="$2"
  export GFORTRAN_UNBUFFERED_ALL=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  track_pid='' monitor_pid=''
  python3 -u visualization/adaptive_dashboard.py --follow "$log" > monitor_events.log 2>&1 & monitor_pid=$!
  # IMPORTANT: tracker accepts a positional path, then --follow boolean.
  python3 -u visualization/watch_second_ring.py "$log" --follow > ring_events.log 2>&1 & track_pid=$!
  cleanup() {
    [[ -z "${monitor_pid:-}" ]] || kill "$monitor_pid" 2>/dev/null || true
    [[ -z "${track_pid:-}" ]] || kill "$track_pid" 2>/dev/null || true
  }
  trap cleanup EXIT
  set +e
  "./$executable" 2>&1 | tee "$log"
  local ret=${PIPESTATUS[0]}
  set -e
  python3 visualization/adaptive_dashboard.py --replay "$log" >/dev/null 2>&1 || true
  cleanup
  trap - EXIT
  echo "Fortran process exit code: $ret (0 does NOT prove end-time convergence)."
  return "$ret"
}
