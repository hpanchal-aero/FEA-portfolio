#!/usr/bin/env bash
# Run CalculiX on <job> inside <dir>: live log, exit code, wall time and peak resident memory.
# Usage: bash scripts/run_ccx.sh simulation/static_v2/medium uCRM9_medium_v2
#        FOLLOW=0 bash scripts/run_ccx.sh ...    (no live log in the terminal)
set -u
CCX="${CCX:-$HOME/miniconda3/envs/ccx223/bin/ccx}"
dir="${1:?usage: run_ccx.sh <dir> <job>}"
job="${2:?usage: run_ccx.sh <dir> <job>}"
cd "$dir" || { echo "ABORT: cannot cd to $dir"; exit 2; }
[ -f "$job.inp" ] || { echo "ABORT: $job.inp not found"; exit 2; }

# line-buffered output so the log is live and survives a crash (C stdio only)
if command -v stdbuf >/dev/null 2>&1; then
  RUN=(stdbuf -oL -eL "$CCX" "$job")
else
  echo "note: stdbuf not found, log will be block-buffered"
  RUN=("$CCX" "$job")
fi

: > ccx.log
start=$(date +%s)
"${RUN[@]}" >> ccx.log 2>&1 &
pid=$!

# peak-memory sampler (VmHWM = high-water mark of resident memory)
(
  peak=0
  while kill -0 "$pid" 2>/dev/null; do
    kb=$(awk '/^VmHWM:/ {print $2}' "/proc/$pid/status" 2>/dev/null)
    if [ -n "$kb" ] && [ "$kb" -gt "$peak" ]; then peak=$kb; fi
    sleep 0.5
  done
  echo "$peak" > .peak_kb
) &
sampler=$!

if [ "${FOLLOW:-1}" = "1" ]; then
  tail -n +1 -f --pid="$pid" ccx.log
fi

wait "$pid"
code=$?
wait "$sampler"
end=$(date +%s)
peak_kb=$(cat .peak_kb 2>/dev/null || echo 0)
rm -f .peak_kb

{
  echo "job                  $job"
  echo "exit code            $code"
  echo "wall time            $((end - start)) s"
  echo "peak resident memory $((peak_kb / 1024)) MiB"
} | tee run_info.txt
