#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
bash visualization/export_all.sh --include-initial --gif "$@"
