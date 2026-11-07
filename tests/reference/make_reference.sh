#!/bin/sh
# make_reference.sh CENTPACKL_DIR [NP] -- regenerate tests/reference/ from CentPackL MHD2d at commit e2626bb.
# The reference is the d_e = 0 contract of IHMHD2d: same input, same rank count, bit-for-bit equal output.
# Input: tests/reference/input = samples/input_gem_1x_de0 (GEM 1x, d_i 1, eta 0.005, eta_H 0, d_e 0,
#        c_whistler 0, t_final 1, dt_out 1). c_whistler = 0 as in every CentPackL reference and paper run.
# Stored: steps.txt (dt of every fluid step and every "EM sub-steps" line of run.log), energy_log.txt,
#         mesh_files/ and the t = 1 frame of every field (*_files/*_p*_001).
# Note: CentPackL MHD2d has no EDRIFT print (that diagnostic was added in the scratch INERTIA2d builds),
#       so the step sequence + energy_log + fields are the reference.
set -e
SRC=${1:?usage: make_reference.sh CENTPACKL_DIR [NP]}
NP=${2:-10}
REF=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d)
git -C "$SRC" archive e2626bb | tar -x -C "$W"
FFTW_DIR=${FFTW_DIR:-/opt/homebrew/opt/fftw}
# MHD2d with the SD2 evolution (the INERTIA2d_SD2 compile line without the INERTIA2d sources)
( cd "$W" && mpicxx -std=c++11 -O3 -Iinclude/ -I$FFTW_DIR/include src/common/disclaimer.cc src/common/minmod3.cc \
    src/common/doublearray*.cc src/2d/common/*.cc src/2d/SD2/*.cc src/2d/samples/MHD2d/*.cc \
    src/2d/sample_mains/main_2d_SD2.cc -o MHD2d )
R=$W/run; mkdir -p "$R"; for s in b1 b2 b3 mesh p rho t u1 u2 u3; do mkdir -p "$R/${s}_files"; done
cp "$REF/input" "$R/input"; cp "$W/MHD2d" "$R/"
( cd "$R" && mpirun -np "$NP" ./MHD2d < /dev/null > run.log 2>&1 )
grep -o 'dt=[^,]*' "$R/run.log" > "$REF/steps.txt" || true
grep 'EM sub-steps' "$R/run.log" >> "$REF/steps.txt" || true
cp "$R/energy_log.txt" "$REF/"
rm -rf "$REF"/*_files; mkdir -p "$REF/mesh_files"; cp "$R"/mesh_files/* "$REF/mesh_files/"
for s in b1 b2 b3 p rho t u1 u2 u3; do mkdir -p "$REF/${s}_files"; cp "$R/${s}_files/"*_001 "$REF/${s}_files/"; done
echo "reference written: $(ls "$REF"/*_files/* | wc -l) files, $(wc -l < "$REF/steps.txt") step lines"
rm -rf "$W"
