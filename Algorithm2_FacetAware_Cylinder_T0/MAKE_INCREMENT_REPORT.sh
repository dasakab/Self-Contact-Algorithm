#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ -f solver_live_windows.log ]]; then f=solver_live_windows.log;else f=solver_live.log;fi
python3 visualization/make_one_line_progress.py "$f" | tee per_increment_report.txt
