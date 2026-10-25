# Local verification of the packaged code (2026-10-25)

- **Linux compilation:** Built current `*.f` + `colsol_sparse.c` with GNU Fortran/SuperLU; build succeeded, with non-fatal legacy Fortran 2018 syntax warnings.
- **Fresh t=0 smoke run:** Started an isolated copy using `bash scripts/run_from_zero_linux.sh --end-time 0.00005 --initial-dt 0.000005`; log verified `ADAPT_START 0 0`, followed by accepted inc1 (45 Newton iterations, t=0.000005, residual below 5e-4) and accepted inc2 (60 iterations, t=0.000010, residual below 5e-4). The smoke run was interrupted in inc3; no long-duration convergence claim is made.
- **Exporter on fresh t=0 smoke outputs:** Produced a reference geometry at t=0 and accepted `.vtu` snapshots inc1 and inc2; PVD index generated successfully.
- **Exporter regression on prior archive:** Generated 10 `.vtu` increments (35–44), PVD series and cross-section GIF from actual archived OUT data plus an **explicitly** selected example checkpoint.
- **Windows MSYS2 UCRT64 native executable:** source and shell scripts prepared; not executable/tested within this Linux environment.
- **Not verified:** full 0.008 or 0.6 second integration, all ring contacts, energy balance and stress accuracy, complete original-vs-current full unit regression.
