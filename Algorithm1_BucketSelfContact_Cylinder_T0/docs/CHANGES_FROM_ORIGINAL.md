# What changed versus the user-supplied `cylinder_250.zip`

Compared **directly** with `cylinder_250.zip/cylinder_250/*.f`, not an inferred historic version.
The full line-by-line exact unified diff is `ORIGINAL_TO_CURRENT.patch`.
Original and packaged SHA-256 hashes are in `source_comparison.json`.

## Fortran files with code differences

- **`read.f`**: initialises `CYLINDER_RAIL=0` before geometry processing, preventing an uninitialised switch from entering the wrong target-mesh branch.
- **`contact2.f`**: substantial experimental contact/self-contact fixes: correct the `VELO` workspace to hold ordinary hitting, inner/outer self-contact hitting **and target nodes** plus target body velocities; avoid target data overwrite; correct `TTF` indexing; repair inner/outer contact-status indexing and sorting; add bounds checks for segment/target searches; resize segment and target-node workspaces; and retain self-contact target lists across increments **inside the same process**. The contact algorithm remains **unvalidated** for physically correct target pairs, signed gaps or normal-force directions.
- **`newrap1.f`**: uses `LOCA/LOCB(NTDOF+1)` for skyline sentinel; sizes the body A scratch geometry for the deformable cylinder rather than the eight target elements; corrects ordinary-contact residual-mask extents; removes the legacy relaxed 0.005-force-residual acceptance; adds adaptive time-step growth/reduction, retry, stricter Newton checks and separately logged inner/outer self-contact. Both inner and outer self-contact storage dimensions are enlarged and persistent target-node arrays are passed to the contact routines.

**Unmodified Fortran files:** `main.f`, `eigensolve.f`, `eulerforint.f`, `force.f`, `modcon.f`, `objstress.f`, `restart.f`, `skyline.f`, `solvers.f`, `stress.f`. Earlier *diagnostic* copies changed `eigensolve.f` and `objstress.f`, but those diagnostic patches are **not included** in this latest solver baseline; do not claim otherwise.

## Non-Fortran additions

- `colsol_sparse.c`: replacement SuperLU sparse direct-solver backend for legacy skyline factorisation (retains the rest of the model code).
- `scripts/build_*`: Linux GNU Fortran and native Windows UCRT64 builds with bounds/FP checks.
- `scripts/run_from_zero_*`: NEW **actual t=0 entry points** (restart option 0 and no checkpoint copy), backup existing results before running.
- `scripts/run_from_checkpoint_*`: deliberately explicit, optional pre-contact diagnostic restart; not the default.
- `visualization/`: independent VTU/PVD/GIF exporter, ParaView cross-section macro, residual monitor and inner/outer ring tracker; fixes the old tracker argument mismatch.
- `UOPT1/OPT_IEXEC`: the sixth line is now `0` (fresh integration, not `RESTART`).
- `UOPT1/O_OEXEC`: default experiment limit time 0.008, maximum 220 iterations, initial dt 2.5e-5; **not** the user's original 0.6-second target.

## Critical limitations and known numerical concerns

1. **Mid-self-contact restart** remains unfaithful: original `RESTART` format doesn't persist all inner/outer contact pair and friction state. Fresh uninterrupted execution is preferred through contact onset.
2. The adaptive retry branch requires independent verification that *all* trial plastic/contact history rolls back upon time-step rejection.
3. No independent proof of force equilibrium, signed contact gap, friction stick/slip accuracy, energy balance, or mesh/time-step convergence for all tube folds.
4. **`YSTRESS_A` state anomaly:** archival inc35 restart uses 1.725e11 Pa yield-state values, inconsistent with the deformable material's initial 3.4e8 Pa yield in `UBIN1/MPROP`; stress results around 1.05 GPa are not yet validated.
5. The fixed Linux continuation reached accepted increment 44, t=0.00310875; no second *physical* fold ring or outer self-contact was validated. Published files are an experimental research code, not a certified impact solver.
