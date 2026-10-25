#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
bash visualization/start_live_monitor.sh "$@"
