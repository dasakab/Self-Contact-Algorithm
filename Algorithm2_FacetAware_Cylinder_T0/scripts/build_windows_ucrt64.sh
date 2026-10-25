#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ "${MSYSTEM:-}" != UCRT64 ]]; then echo 'Run from MSYS2 UCRT64, not normal MSYS, CMD or PowerShell.' >&2; exit 2; fi
for cmd in gcc gfortran python3;do command -v "$cmd" >/dev/null || { echo "Missing $cmd" >&2; exit 3; };done
[[ -f /ucrt64/include/superlu/slu_ddefs.h ]] || { echo 'Install mingw-w64-ucrt-x86_64-superlu' >&2; exit 4; }
mkdir -p build/windows
flags=(-O2 -g -fcheck=bounds -fbacktrace -ffpe-trap=invalid,zero,overflow -ffixed-line-length-none -fallow-argument-mismatch)
files=(main read contact2 eigensolve eulerforint force modcon newrap1 objstress restart skyline stress facet_search)
objs=()
for f in "${files[@]}";do
  echo "Compiling $f.f"
  gfortran -c "${flags[@]}" "$f.f" -o "build/windows/$f.o"
  objs+=("build/windows/$f.o")
done
grep -Eq '^ *SUBROUTINE MATIN\(' solvers.f || { echo 'Expected MATIN boundary not found in solvers.f'; exit 5; }
sed -n '/^ *SUBROUTINE MATIN(/,$p' solvers.f > build/windows/solvers_nocol.f
gfortran -c "${flags[@]}" build/windows/solvers_nocol.f -o build/windows/solvers_nocol.o
gcc -O2 -c colsol_sparse.c -o build/windows/colsol_sparse.o
echo 'Linking solver_windows.exe with SuperLU'
gfortran -g -O2 -o solver_windows.exe "${objs[@]}" build/windows/solvers_nocol.o build/windows/colsol_sparse.o -lsuperlu -lm
echo 'Built native Windows executable: solver_windows.exe'
