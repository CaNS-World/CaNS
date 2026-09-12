#!/usr/bin/env bash
set -euo pipefail
testdir=$(cd "$(dirname "$0")" && pwd)
rundir="$testdir/../../run"
mkdir -p "$rundir"
python3 "$testdir/test.py" 2>&1 | tee "$rundir/tgv.log"
