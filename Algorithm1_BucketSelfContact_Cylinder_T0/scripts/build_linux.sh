#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
for cmd in gcc gfortran;do command -v "$cmd" >/dev/null || { echo "Missing $cmd" >&2; exit 3; };done
mkdir -p build/linux
flags=(-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=invalid,zero,overflow -ffixed-line-length-none -fallow-argument-mismatch)
files=(main read contact2 eigensolve eulerforint force modcon newrap1 objstress restart skyline stress)
objs=()
for f in "${files[@]}";do
  echo "Compiling $f.f"
  gfortran -c "${flags[@]}" "$f.f" -o "build/linux/$f.o"
  objs+=("build/linux/$f.o")
done
grep -Eq '^ *SUBROUTINE MATIN\(' solvers.f || { echo 'Expected MATIN boundary not found in solvers.f'; exit 5; }
sed -n '/^ *SUBROUTINE MATIN(/,$p' solvers.f > build/linux/solvers_nocol.f
gfortran -c "${flags[@]}" build/linux/solvers_nocol.f -o build/linux/solvers_nocol.o
gcc -O2 -c colsol_sparse.c -o build/linux/colsol_sparse.o
echo 'Linking solver_linux with SuperLU'
gfortran -g -O2 -o solver_linux "${objs[@]}" build/linux/solvers_nocol.o build/linux/colsol_sparse.o -lsuperlu -lm
echo 'Built Linux executable: solver_linux'
