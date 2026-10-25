# Algorithm2 — Facet-aware self-contact for a deforming cylinder

**Status: experimental research code, NOT validated nonpenetrating tube collapse.**

This is the complete source/batch/post-processing folder for the **Algorithm2**
facet-aware trial, with a **genuine fresh start at t=0** by default.
The underlying structural/contact solver is Fortran with a SuperLU sparse linear
solver written in C. `facet_search.f` augments the earlier contact candidate
search and closest-point projection. The legacy Lagrange-multiplier contact
active set is still present; a guarantee of zero interpenetration or detection
of all folds has **not** been established.

## Quick start — Windows 11, MSYS2 UCRT64

1. Unzip the whole package into a **NEW** folder, e.g.
   `C:\Users\Admin\Downloads\thesis\Algorithm2_FacetAware_Cylinder_T0`.
   Do not merge it into the old `cylinder_250_second_ring` folder.
2. Open **MSYS2 UCRT64** (the prompt must display `UCRT64`).
3. If needed, install the compiler, C compiler, SuperLU and Python:

   ```bash
   pacman -S --needed mingw-w64-ucrt-x86_64-gcc \
     mingw-w64-ucrt-x86_64-gcc-fortran \
     mingw-w64-ucrt-x86_64-superlu \
     mingw-w64-ucrt-x86_64-python
   ```
4. Navigate to the **actual extracted folder** (adjust location as needed):

   ```bash
   cd /c/Users/Admin/Downloads/thesis/Algorithm2_FacetAware_Cylinder_T0
   pwd
   ls facet_search.f UOPT1/OPT_IEXEC scripts/run_from_zero_ucrt64.sh
   ```
5. Compile in the first terminal:

   ```bash
   bash scripts/build_windows_ucrt64.sh
   ```

   Success ends with `Built native Windows executable: solver_windows.exe`.
6. Start **from time ZERO**:

   ```bash
   bash scripts/run_from_zero_ucrt64.sh
   ```

   This sets record 6 of `UOPT1/OPT_IEXEC` to **0**, archives any old OUT/RESTART
   under `run_backups/`, and starts the solver without importing increment 35.
   Check for `ADAPT_START` reporting increment 0 and time 0, then increment 1.
   You can also build-and-run in one command: `bash START_WINDOWS_UCRT64.sh`.

   Optional: `bash scripts/run_from_zero_ucrt64.sh --initial-dt 5e-6 --max-iterations 220 --end-time 0.008`

## Live convergence dashboard (second terminal)

```bash
cd /c/Users/Admin/Downloads/thesis/Algorithm2_FacetAware_Cylinder_T0
bash START_LIVE_DASHBOARD.sh
```

Open `http://localhost:8765/monitor.html` on the **same PC**. This reads the
numerical log on **every** refresh, avoiding stale HTTP 304 snapshots. The
page shows Newton residuals, tolerance, contact counts, current increment and
last accepted increment. If the solver terminates, it cannot show additional
numeric progress. If the port is already used, stop the old local server first.

Raw numerical output: `tail -f solver_live_windows.log` in another terminal.
Stop log-following with Ctrl+C; that does not stop a solver in a different terminal.

## ParaView time-series and cross-section GIF

Install visualisation dependencies in MSYS2 UCRT64:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-python-matplotlib \
  mingw-w64-ucrt-x86_64-python-pillow
```

After the solver has written **complete** steps to `OUT`, run:

```bash
bash EXPORT_PARAVIEW_GIF.sh
```

Or run VTU/PVD export without the GIF dependencies:

```bash
bash visualization/export_all.sh --include-initial
```

The results go into `PARAVIEW_EXPORT/`:

- `cylinder_inc0000.vtu`: **reference** shape (not a computed increment),
  placed with the same model origin as the solver.
- `cylinder_incXXXX.vtu`: one dataset per **complete saved increment**.
- `cylinder_timesteps.pvd`: ParaView timeline (open this to animate).
- `cylinder_cross_section.gif`: animation using a **fixed X-Z slice at y=0**.
- `PNG_frames/frame_XXXX.png`: separate image for each saved increment.

Missing/failed increments are not invented or interpolated. Export only after
results have been saved; rerun the exporter later to include more increments.
The GIF uses a fixed 0–1100 MPa legend range for comparability; it does not
imply the abnormal stress values have been validated.

### Inside ParaView

1. File -> Open `PARAVIEW_EXPORT/cylinder_timesteps.pvd`; click **Apply**.
2. Use **Filters -> Slice**, with **Origin** `(0,0,22)` and **Normal** `(0,1,0)`.
3. Set **Surface With Edges**. Set camera to view along Y, Z vertical.
4. Colour using `equivalent_stress_Pa`, `displacement_magnitude`, or
   `equivalent_plastic_strain`.
5. Click **Play** and use **File -> Save Animation** for another video/GIF.
6. Optional: execute `visualization/paraview_cross_section_macro.py` in ParaView
   **View -> Python Shell -> Run Script** after loading the `.pvd`.
   The Python shell is inside ParaView; do not execute this macro with standard Python.

## Surface penetration audit (optional Python dependencies)

Install VTK, NumPy, SciPy using a compatible Python environment (e.g. Linux/WSL
venv) and run:

```bash
python3 visualization/penetration_audit.py PARAVIEW_EXPORT/cylinder_inc0037.vtu
```

Audit counts are **triangle-pair intersections**, not distinct physical contacts.
The legacy contact active-set and friction handling still require physical checks.

## One-sentence accepted increment report

```bash
bash MAKE_INCREMENT_REPORT.sh
```

It reads the solver log and writes `per_increment_report.txt` with accepted
Newton increments. The live dashboard uses the same underlying log.

## Linux / WSL2

Install `gfortran`, `gcc`, `libsuperlu-dev`, and Python 3. Then run:

```bash
bash scripts/build_linux.sh
bash scripts/run_from_zero_linux.sh
bash visualization/export_all.sh --include-initial --gif
```

GNU Fortran compiler/linker warnings about Fortran 2018 legacy syntax are not
by themselves errors. However, floating-point exceptions, bounds failures,
NaN residuals and geometric crossings mean a failed or unvalidated simulation.

## Numerical defaults

- Body A: 2,640 8-node hexahedral elements, 4,644 nodes.
- Body B: 8-element rigid/obstacle mesh in project inputs.
- Impact initial velocity from `UOPT1/DYNA_IN` (Z = -250 in model units).
- Initial Δt = `5e-6`; iteration limit = `220`; convergence tolerance ≈ `5e-4`.
- Nominal analysis end-time = `0.008`. Later numerical collapse unverified.
- Force tolerance is not the same as a no-penetration guarantee.

## Contact changes and attribution

- `contact2.f`: calls experimental facet candidates on inner and outer surfaces,
  and uses closest-point projection in segment selection.
- `facet_search.f`: searches expanded facet bounding boxes, computes clamped
  quadrilateral closest-point projection, checks local surface orientation and
  filters candidates using signed gaps.
- Other earlier source fixes include skyline bounds, residual masking, node
  indexing/VELOCITY state, adaptive time stepping, and sparse factorisation.
- The `BUCKET` routine in the original source is a cylindrical Z–theta ordering,
  **not** a generic 3D Cartesian spatial hash. The educational posters are
  conceptual diagrams, not exact computational schematics.
- The structural/contact solver still uses Lagrange multipliers and legacy
  active-set logic. The illustrative penalty-contact equation in some posters
  is **not** proof this code uses a penalty law.

See `docs/README_FACET_CONTACT.md` and the source changes in `docs/*.patch`.
The pre-contact increment-35 checkpoint in `examples/checkpoints` is only for
independent comparison; **fresh-start scripts do not load it**.

## Test status

This folder is assembled from the latest compiled Algorithm2 trial sources,
not from an independently validated final-collapse solution. The earlier
Linux check from the increment-35 checkpoint converged increments 36 and 37
with no boundary triangle crossings reported by the audit, then stopped
without accepting increment 38. A full **Algorithm2** run starting at t=0 has
not been demonstrated through the second self-contact ring.

Please retain credit for the original self-contact author Sarma, and verify
permissions before publishing source code under an open-source licence.
