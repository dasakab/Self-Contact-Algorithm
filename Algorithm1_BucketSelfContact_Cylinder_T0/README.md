# Algorithm1 — Bucket / cylindrical-grid self-contact (cylinder impact)

**Research software; not validated to prevent self-penetration.**

This is the reproducible **Algorithm1** counterpart to `Algorithm2_FacetAware_Cylinder_T0`.
It is packaged with the same mesh, material inputs, nonlinear solver, sparse SuperLU
backend, adaptive time stepping, fresh-start scripts, VTU/PVD/GIF exporter and live
monitor. The key intentional difference is `contact2.f`: **Algorithm1 retains
Sarma's cylindrical Z–theta `BUCKET` ordering, `SNODE_HT_IN/OUT`, legacy nearest-node/
segment contact candidate search and contact assembly**. It does not call
`SELF_FACET_CANDIDATES`, `SAFE_QUAD_PROJECT`, or `facet_search.f`.

The `contact2.f` file is *not byte-identical* to the supplied original. Earlier
memory-safety, contact array workspace and persistence fixes have been applied.
See `docs/CHANGES_FROM_ORIGINAL.md` and the source diff
`docs/ORIGINAL_TO_CURRENT.patch` for the differences.

## 1 — Native Windows / MSYS2 UCRT64 (recommended for your PC)

Extract this ZIP **to its own directory**; do not merge with the Algorithm2 files.
Open MSYS2 UCRT64 (the prompt must show `UCRT64`).

Install dependencies if not already installed:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-gcc \
  mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-superlu \
  mingw-w64-ucrt-x86_64-python
```

Example location (change if you extracted elsewhere):

```bash
cd /c/Users/Admin/Downloads/thesis/Algorithm1_BucketSelfContact_Cylinder_T0
bash scripts/build_windows_ucrt64.sh
bash scripts/run_from_zero_ucrt64.sh
```

Alternatively, build + start with `bash START_WINDOWS_UCRT64.sh`.
A successful build creates `solver_windows.exe` in the project root. A true
fresh run shows `ADAPT_START 0 0 ...` before starting increment 1. The script
sets record 6 of `UOPT1/OPT_IEXEC` to `0`, archives previous `RESTART` and `OUT`
in `run_backups/`, and starts without loading a checkpoint.

To set run parameters, use e.g.

```bash
bash scripts/run_from_zero_ucrt64.sh --initial-dt 5e-6 --max-iterations 220 --end-time 0.008
```

Do not directly invoke the old checkpoint script if you want a fresh start.
For full historical 0.6-second nominal end time, use `--end-time 0.6`;
reaching that time has NOT been verified.

## 2 — Linux / WSL2

```bash
sudo apt-get update
sudo apt-get install gfortran gcc libsuperlu-dev python3
bash scripts/build_linux.sh
bash scripts/run_from_zero_linux.sh
```

## 3 — Live Newton convergence dashboard

In a second MSYS2 UCRT64 (or Linux/WSL2) terminal, from the project folder:

```bash
bash START_LIVE_DASHBOARD.sh
```

Visit `http://localhost:8765/monitor.html` locally. The monitor reads fresh
solver logs on request. Raw output can be watched in another terminal:

```bash
tail -f solver_live_windows.log
```

The dashboard does not start or resume the solver. If there is no new log
because the solver has stopped, it will show the last recorded state.

## 4 — ParaView and GIF post-processing

On Windows MSYS2 UCRT64, install optional plotting packages:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-python-matplotlib \
  mingw-w64-ucrt-x86_64-python-pillow
```

Run after complete converged increments have been written to `OUT`:

```bash
bash EXPORT_PARAVIEW_GIF.sh
```

For VTU/PVD export only (no matplotlib/Pillow required):

```bash
bash visualization/export_all.sh --include-initial
```

Files appear in `PARAVIEW_EXPORT/`:

- `cylinder_inc0000.vtu` — **reference** at t=0 (not a computed increment)
- `cylinder_incXXXX.vtu` — available saved deformed increments
- `cylinder_timesteps.pvd` — ParaView animation index
- `cylinder_cross_section.gif` — X–Z cross section at y=0, fixed camera
- `PNG_frames/` — per-increment animation frames

In ParaView, open `cylinder_timesteps.pvd`, Apply, Slice with normal `(0,1,0)`
and origin `(0,0,22)`, Surface With Edges, color using
`equivalent_stress_Pa` or `equivalent_plastic_strain`, then Play. The
included `visualization/paraview_cross_section_macro.py` can set the camera
from ParaView's Python Shell (not the ordinary terminal).

The plotting pipeline includes a solver-reference-coordinate offset for the
undeformed t=0 frame. Missing increments are never interpolated.

## 5 — Convergence log & geometry audit

```bash
bash MAKE_INCREMENT_REPORT.sh
python3 visualization/penetration_audit.py PARAVIEW_EXPORT/cylinder_inc0037.vtu
```

The second command requires compatible `vtk`, NumPy and SciPy libraries.
The audit reports nonincident surface-triangle intersections, not a validated
measurement of physical gap or independent self-contact pairs. More detailed
signed gaps, active multipliers, mesh and energy checks remain necessary.

## 6 — Matching Algorithm2 settings

For a controlled A/B experiment use the exact same input geometry, mesh,
material settings, initial velocity, end time, Δt, convergence tolerance and
compiler settings. This package uses:

| Parameter | Packaged Algorithm1 value |
|---|---|
| Start time | t=0 |
| Restart flag | 0 |
| Initial Δt | 5e-6 s |
| Max Newton iterations | 220 |
| Force residual tolerance | approx. 5e-4 |
| Nominal end time | 0.008 s |

Both algorithms contain adaptive time-step rules; their *accepted time grids*
may diverge after contact onset. Compare their outputs at matched **physical
simulation times**, not simply increment numbers. Any acceptance must also
pass independent surface-intersection, energy and mesh checks.

## 7 — What is and is not preserved from Sarma's source

This experimental baseline preserves the legacy Z–theta node ordering and
node/segment-based self-contact search. Earlier pragmatic corrections to
`contact2.f` and `newrap1.f` are retained so it compiles and progresses.
`colsol_sparse.c` adds a SuperLU sparse solver. The original `BUCKET`
method is **not** a general-purpose Cartesian 3D spatial grid/hash; educational
illustrations of spatial cubes are a *conceptual analogy*, not a literal
representation of this implementation.

Prior comparison snapshots from increments 35–44, if provided in `examples/`,
come from an earlier experimental adaptive run, not this exact matched-Δt
configuration and not a validated contact solution. Do not represent those
as independent Algorithm1-vs-Algorithm2 validation.

## 8 — Warnings and author credit

The legacy Lagrange-multiplier active-set/self-contact logic can still miss
facet crossings and produce penetrations. There is no claim of guaranteed
zero penetration, correct second-ring contact, or validated absolute stress.
The archived increment-35 `YSTRESS_A` field has a material-state discrepancy;
previous approximately 1.05 GPa stress peaks need independent investigation.

The original self-contact implementation is attributed by the project owner
to **Sarma**. Original source comments did not provide a formal author header.
Keep attribution and obtain all necessary redistribution/licensing permissions
before publishing the complete source to GitHub. `docs/LICENSING.md` documents
that review; no open-source licence is granted by this archive.
