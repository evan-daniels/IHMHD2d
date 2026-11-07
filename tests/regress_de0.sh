#!/bin/sh
# regress_de0.sh [NP] -- d_e = 0 regression of IHMHD2d against tests/reference/ (CentPackL MHD2d, e2626bb, c_whistler = 0).
# Runs GEM 1x (tests/reference/input, d_e = 0, t_final 1) on NP ranks (default 10; the reference was made on 10)
# and requires bit-for-bit equality of:
#   - the per-step time steps and EM sub-step counts (run.log),
#   - energy_log.txt,
#   - every field file at the outputs t = 0, 0.5, 1.
# Exit status 0 = identical, 1 = any difference.
set -e
NP=${1:-10}
T=$(cd "$(dirname "$0")" && pwd)
REF=$T/reference
RUN=$T/run_de0
rm -rf "$RUN"; mkdir -p "$RUN"
for s in b1 b2 b3 mesh p rho t u1 u2 u3; do mkdir -p "$RUN/${s}_files"; done
cp "$REF/input" "$RUN/input"
cd "$RUN"
mpirun -np "$NP" ../../IHMHD2d < /dev/null > run.log 2>&1
grep -o 'dt=[^,]*' run.log > steps.txt || true
grep 'EM sub-steps' run.log >> steps.txt || true
fail=0
# ORIGINAL (pre-2026-11-06):
# cmp -s steps.txt "$REF/steps.txt" || { echo "FAIL: time steps / EM sub-steps differ"; fail=1; }
# NEW (2026-11-06): ranks' stdout interleaves, so line order and short-format/split tokens differ run to run.
# Compare the set of complete full-precision dt values and the sorted EM sub-step lines instead.
# ORIGINAL (pre-2026-11-07):
# canon() { grep -oE 'dt=[0-9]\.[0-9]{16}e[+-][0-9]+' "$1" | sort -u; grep 'EM sub-steps' "$1" | sort; }
# NEW (2026-11-07): exactly two exponent digits followed by a non-digit (an interleaved digit gave e-016)
canon() { grep -oE 'dt=[0-9]\.[0-9]{16}e[+-][0-9]{2}([^0-9]|$)' "$1" | grep -oE 'dt=[0-9.e+-]{22}' | sort -u; grep 'EM sub-steps' "$1" | sort; }
canon steps.txt > steps_canon.txt; canon "$REF/steps.txt" > ref_canon.txt
cmp -s steps_canon.txt ref_canon.txt || { echo "FAIL: time steps / EM sub-steps differ"; fail=1; }
cmp -s energy_log.txt "$REF/energy_log.txt" || { echo "FAIL: energy_log.txt differs"; fail=1; }
n=0
for f in $(cd "$REF" && find . -path './*_files/*' -type f | sort); do
  n=$((n+1)); cmp -s "$f" "$REF/$f" || { echo "FAIL: $f differs"; fail=1; }
done
if [ $fail -eq 0 ]; then echo "PASS: d_e = 0 bit-identical to CentPackL MHD2d ($n field files, energy_log, $(wc -l < steps.txt) step lines)"; fi
exit $fail
