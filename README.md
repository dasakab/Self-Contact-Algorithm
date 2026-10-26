# Self-Contact-Algorithm

**Nonlinear finite-element cylinder impact and experimental self-contact solvers**

This repository contains two related Fortran finite-element solver variants for deforming-cylinder impact calculations. Both use a C/SuperLU sparse linear-solver backend, Newton iterations, adaptive time stepping, and scripts for run monitoring and ParaView export.



## 1. Algorithms

| | Algorithm1 | Algorithm2 |
|---|---|---|
| Directory | `Algorithm1_BucketSelfContact_Cylinder_T0/` | `Algorithm2_FacetAware_Cylinder_T0/` |
| Self-contact search | Legacy cylindrical Z–theta `BUCKET` ordering and node/segment candidate search | Experimental facet-aware candidates and closest-point projection |
| Main distinguishing code | `contact2.f` | `contact2.f` plus `facet_search.f` |
| Linear algebra | Fortran + C/SuperLU | Fortran + C/SuperLU |
| Fresh-start and monitoring scripts | Included | Included |

Algorithm1 retains the original *style* of Sarma's contact search, with subsequent bug fixes; it is not byte-identical to the original source. Algorithm2 modifies contact candidate selection, but retains aspects of the legacy Lagrange-multiplier active-set implementation. Neither guarantees detection of all self-contact events.

Keep the folders separate; each contains its own input, build, `OUT/`, `RESTART/`, and log files.

## 2. Download

Open [the GitHub repository](https://github.com/dasakab/Self-Contact-Algorithm), select **Code → Download ZIP**, extract it, and keep the two solver folders intact. Git is optional for downloading; a compiler and dependencies are required to **build and run** the solvers.

Expected layout (selected files shown):

```text
Self-Contact-Algorithm/
├── README.md
├── Algorithm1_BucketSelfContact_Cylinder_T0/
│   ├── contact2.f
│   ├── colsol_sparse.c
│   ├── UOPT1/
│   ├── scripts/
│   ├── visualization/
│   ├── START_WINDOWS_UCRT64.sh
│   ├── START_LINUX_WSL.sh
│   ├── START_LIVE_DASHBOARD.sh
│   ├── MAKE_INCREMENT_REPORT.sh
│   └── EXPORT_PARAVIEW_GIF.sh
└── Algorithm2_FacetAware_Cylinder_T0/
    ├── contact2.f
    ├── facet_search.f
    ├── colsol_sparse.c
    ├── UOPT1/
    ├── scripts/
    ├── visualization/
    ├── START_WINDOWS_UCRT64.sh
    ├── START_LINUX_WSL.sh
    ├── START_LIVE_DASHBOARD.sh
    ├── MAKE_INCREMENT_REPORT.sh
    └── EXPORT_PARAVIEW_GIF.sh
```

## 3. Windows 10/11: install dependencies

The supplied native-Windows scripts require **MSYS2 UCRT64** (not PowerShell, Command Prompt, or the plain MSYS shell).

1. Install MSYS2 from [msys2.org](https://www.msys2.org/).
2. Open **MSYS2 UCRT64** from the Start menu. Confirm the prompt shows `UCRT64`.
3. Update the package database and packages:

   ```bash
   pacman -Syu
   ```

   If MSYS2 requests that you close the terminal, reopen **MSYS2 UCRT64** and run `pacman -Syu` again.
4. Install the required compilers, sparse solver, and Python:

   ```bash
   pacman -S --needed \
     mingw-w64-ucrt-x86_64-gcc \
     mingw-w64-ucrt-x86_64-gcc-fortran \
     mingw-w64-ucrt-x86_64-superlu \
     mingw-w64-ucrt-x86_64-python
   ```

5. Check the tools:

   ```bash
   echo "$MSYSTEM"             # UCRT64
   gfortran --version
   gcc --version
   python3 --version
   ```

Paths such as `C:\Users\USERNAME\Downloads` are written as `/c/Users/USERNAME/Downloads` in MSYS2. Replace the example paths below with your actual extracted folder.

### 3.1 Build and run Algorithm1 from t = 0

In an **MSYS2 UCRT64** terminal:

```bash
cd /c/Users/USERNAME/Downloads/Self-Contact-Algorithm-main/Algorithm1_BucketSelfContact_Cylinder_T0
bash scripts/build_windows_ucrt64.sh
bash scripts/run_from_zero_ucrt64.sh \
  --initial-dt 5e-6 --max-iterations 220 --end-time 0.008
```

Expected compiler output ends with `Built native Windows executable: solver_windows.exe`. The fresh-start script should identify an initial increment/time of zero and print an `ADAPT_START` marker.

**Alternative:** `bash START_WINDOWS_UCRT64.sh` builds and starts with the configured input defaults. Do not run both launch methods for the same job.

### 3.2 Build and run Algorithm2 from t = 0

Open a **separate MSYS2 UCRT64** terminal:

```bash
cd /c/Users/USERNAME/Downloads/Self-Contact-Algorithm-main/Algorithm2_FacetAware_Cylinder_T0
bash scripts/build_windows_ucrt64.sh
bash scripts/run_from_zero_ucrt64.sh \
  --initial-dt 5e-6 --max-iterations 220 --end-time 0.008
```

The expected executable is also `solver_windows.exe`, but **inside the Algorithm2 folder**. Its build additionally compiles `facet_search.f`.

**Alternative:** `bash START_WINDOWS_UCRT64.sh`.

> The packaged Windows build scripts have been inspected; successful execution on every Windows configuration or completion of the full impact simulation has **not** been verified.

## 4. Linux / Ubuntu / WSL2

On Ubuntu or Ubuntu under WSL2, install:

```bash
sudo apt-get update
sudo apt-get install -y gfortran gcc libsuperlu-dev python3
```

From the extracted repository root, use **one solver folder at a time**:

**Algorithm1**

```bash
cd Algorithm1_BucketSelfContact_Cylinder_T0
bash scripts/build_linux.sh
bash scripts/run_from_zero_linux.sh \
  --initial-dt 5e-6 --max-iterations 220 --end-time 0.008
```

**Algorithm2** (from the repository root in another terminal):

```bash
cd Algorithm2_FacetAware_Cylinder_T0
bash scripts/build_linux.sh
bash scripts/run_from_zero_linux.sh \
  --initial-dt 5e-6 --max-iterations 220 --end-time 0.008
```

Both scripts create `solver_linux` in the **respective** project directory. `bash START_LINUX_WSL.sh` is a build-and-fresh-run shortcut.

## 5. Configure and start a job

The recommended launch scripts configure `UOPT1/OPT_IEXEC` and `UOPT1/O_OEXEC` before executing the solver.

| Option | Typical packaged value | Meaning |
|---|---|---|
| `--initial-dt` | `5e-6` | Initial physical time increment |
| `--max-iterations` | `220` | Maximum Newton iterations |
| `--end-time` | `0.008` | Requested simulation end time |
| Force residual tolerance | Approximately `5e-4` | Solver convergence setting |

Example:

```bash
bash scripts/run_from_zero_ucrt64.sh \
  --initial-dt 5e-6 --max-iterations 220 --end-time 0.008
```

For Linux substitute `run_from_zero_linux.sh`.

The **fresh-start** script sets the restart flag to `0` (record 6 of `UOPT1/OPT_IEXEC`), moves earlier results and checkpoints into `run_backups/`, then starts with **no previous checkpoint loaded**. The nominal end time is a target, **not** a claim that the calculation will reach it. Adaptive stepping means the two algorithms may accept different time increments.

Check the original mesh, material, and loading files under `UOPT1/` before changing the physical problem. The example cylinder package has 2,640 hexahedral elements on the deformable body, plus an obstacle mesh; the packaged initial impact velocity is configured in `UOPT1/DYNA_IN`.

## 6. Monitor a running job

Keep the solver running in **Terminal 1**. Open **Terminal 2** (same operating environment), change to the **same algorithm folder**, and start:

```bash
bash START_LIVE_DASHBOARD.sh
```

Open **http://localhost:8765/monitor.html** in your browser, on the same computer. The local dashboard reads the solver log and shows the latest recorded Newton convergence, iterations, accepted increment, and contact-related diagnostics.

For a live text log in **Terminal 3**:

**Windows / MSYS2 UCRT64**

```bash
tail -f solver_live_windows.log
```

**Linux / WSL2**

```bash
tail -f solver_live.log
```

Stop `tail -f` using `Ctrl+C`; this stops **log viewing**, not the solver in Terminal 1.

### How to interpret progress

- `ADAPT_START`: starting increment, time, and adaptive-step state. Check that a fresh job starts at zero.
- `INC NO` / `ITR NO`: increment being attempted and Newton iteration.
- `UNBALANCE FORCE`: iterative force-residual diagnostic.
- `ADAPT_ACCEPT`: solver reports an accepted increment; use its physical time when comparing runs.
- Contact candidate counts or ring markers: diagnostic detections, **not proof of a separate physical contact ring**.
- Missing new log entries: the program may have completed, stopped, crashed, or be stalled. Check Terminal 1 and the executable's exit message.

You can also generate a readable progress summary:

```bash
bash MAKE_INCREMENT_REPORT.sh
```

This writes `per_increment_report.txt`. A pending increment is **not** a converged result.

**Two dashboards:** both solvers default to port `8765`. To monitor them simultaneously on one PC, use separate isolated sessions or modify one dashboard's port configuration. Otherwise start **only one local dashboard at a time**, and change the current working directory when switching algorithms.

## 7. Results and ParaView visualisation

The solver stores native increment data in `OUT/`, and restart data in `RESTART/` when generated. To convert **complete, saved increments** into ParaView files, run this command from **either algorithm folder**:

```bash
bash visualization/export_all.sh --include-initial
```

For an exported time series **plus GIF**, install the optional packages in MSYS2 UCRT64:

```bash
pacman -S --needed \
  mingw-w64-ucrt-x86_64-python-matplotlib \
  mingw-w64-ucrt-x86_64-python-pillow
```

Then run:

```bash
bash EXPORT_PARAVIEW_GIF.sh
```

The generated `PARAVIEW_EXPORT/` folder can contain:

| File | Purpose |
|---|---|
| `cylinder_inc0000.vtu` | Undeformed **reference** configuration; not a calculated increment |
| `cylinder_incXXXX.vtu` | Exported saved increment |
| `cylinder_timesteps.pvd` | ParaView time-series index |
| `cylinder_cross_section.gif` | Fixed X–Z cross-section animation |
| `PNG_frames/` | Per-increment animation frames |

**ParaView steps**

1. Open `PARAVIEW_EXPORT/cylinder_timesteps.pvd` and click **Apply**.
2. Colour by `equivalent_stress_Pa`, `equivalent_plastic_strain`, or `displacement_magnitude` as available.
3. For a mid-section, apply a **Slice** with Normal `(0, 1, 0)` and Origin `(0, 0, 22)` for the packaged example geometry.
4. Select **Surface With Edges** if you want to inspect mesh deformation.
5. Use the animation controls to inspect saved time steps.

The macro `visualization/paraview_cross_section_macro.py` runs from **ParaView's Python Shell**, not ordinary Python. The exporter does not invent missing or failed increments.

## 8. Self-contact and convergence checks

Run the one-line report:

```bash
bash MAKE_INCREMENT_REPORT.sh
```

For an independent *surface-triangle intersection* diagnostic, first ensure `vtk`, `numpy`, and `scipy` are installed in the **Python environment used for this command**. Then run, replacing `XXXX` with an actual saved increment:

```bash
python3 visualization/penetration_audit.py \
  PARAVIEW_EXPORT/cylinder_incXXXX.vtu
```

This is an **audit**, not complete contact validation. Triangle-pair intersections do not equal independent physical contact pairs, and zero reported intersections alone does not prove correct contact forces.

For an Algorithm1/Algorithm2 comparison, match the geometry, mesh, material, initial velocity, convergence settings and compiler options; compare the solutions at matching **physical times**, not necessarily matching increment numbers. Inspect force residuals, energy balance, mesh quality, penetration, and stress plausibility before treating the results as physically trustworthy.

## 9. Restarting and protecting results

Use `run_from_zero_*` only when you intentionally want to start over. It archives prior `OUT/`, `RESTART/`, logs, and visualisation data into a timestamped directory under `run_backups/`.

Both packages also provide experimental checkpoint launchers, e.g.:

```bash
# Algorithm1 on Windows
bash scripts/run_from_checkpoint_ucrt64.sh examples/checkpoints/INC35_VERIFIED

# Algorithm2 on Windows
bash scripts/run_from_checkpoint_ucrt64.sh examples/checkpoints/INC35_PRECONTACT
```

For Linux use `scripts/run_from_checkpoint_linux.sh` with the corresponding checkpoint path.

**Important:** The legacy checkpoint format does **not** fully save the active self-contact pair set. Restarting after contact begins may not reproduce an uninterrupted run. Do not claim exact restart equivalence or delete valuable checkpoints. For an uninterrupted run, keep the solver process running through the contact event.

## 10. Troubleshooting

| Symptom | Action |
|---|---|
| `Run from MSYS2 UCRT64` | Open the **MSYS2 UCRT64** terminal, not plain MSYS or PowerShell. |
| `gfortran: command not found` | Install the UCRT64 GCC/GFortran packages; check `echo "$MSYSTEM"` and `gfortran --version`. |
| Missing `slu_ddefs.h` or `-lsuperlu` | Install the UCRT64 SuperLU package on Windows, or `libsuperlu-dev` on Ubuntu. |
| `scripts/...: No such file` | `cd` into the **individual algorithm directory**, not the repository root. |
| Dashboard says address already in use | Stop the previous monitor that occupies port `8765`; the other solver may still be running. |
| Dashboard does not advance | Check the underlying solver log, Terminal 1, and whether any further accepted steps were written. |
| Bounds failure, floating-point exception, NaN, or segmentation fault | Treat as a solver failure; preserve logs and checkpoints for debugging. Do not assume the last attempted increment converged. |
| No VTU files | Check `OUT/` for completed saved increments; rerun the export script. |
| GitHub repository too large | Exclude generated executables, logs, `OUT/`, `RESTART/`, `run_backups/`, and large `PARAVIEW_EXPORT/` output from ordinary source-code commits. |

## 11. Known limitations and test status

- This project has **not** demonstrated validated full-cylinder collapse or guaranteed elimination of self-penetration.
- Algorithm1's BUCKET search uses a cylindrical grid/order; it is **not** a general Cartesian 3D spatial hash.
- Algorithm2 adds facet-based search and projection but retains experimental search restrictions and legacy contact active-set logic.
- The iterative solver can converge while the exported surfaces still intersect.
- The historic checkpoint-35/step-36–37 trials are not proof that a newly compiled fresh-start job will converge to the same late time state.
- Reported large stresses, including historical values near 1 GPa, remain subject to material-state and numerical verification.
- Source scripts and documentation have been examined; full native-Windows execution and end-time convergence have **not** been independently verified for both current packages.

## 12. Attribution and licensing

The underlying self-contact implementation is attributed by the project owner to **Sarma**. The repository contains experimental modifications and debugging changes. Keep the original attribution and consult the `docs/` files inside each solver folder for modification history.

**No open-source licence is implied merely because the repository is public.** Confirm ownership and redistribution rights before publishing or accepting contributions.

---

**Suggested workflow:** Build Algorithm1 → run from zero → monitor and export → repeat independently for Algorithm2 → compare **accepted physical times**, intersection audits, and numerical convergence.
