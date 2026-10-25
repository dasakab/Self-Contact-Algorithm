# Cylinder impact: experimental facet-contact correction (25 October 2026)

## Status — NOT a full nonpenetrating solution

This is an independently compiled diagnostic branch; **the requested complete fold collapse, second independent ring and zero-penetration guarantee have NOT been achieved**. The solver still uses the original Lagrange-multiplier/force heuristic in `CHECK_CONTACT`. It is not a production contact implementation.

## Verified tests, from original pre-contact increment-35 checkpoint

Run A (`DT=1.875e-5`):
- Increment 36 converged, no self-contact candidates accepted.
- Increment 37 converged at t=0.0030875, but an independent boundary triangle-triangle audit reported **192 nonincident triangle intersections**. This is a failed geometric verification despite Newton convergence.
- Increment 38 identified 48 penetrating inner node-to-facet candidates but was not accepted; the bounded run ended.

Run B (`DT=5e-6`, packaged setting):
- Increment 35: saved original precontact starting checkpoint, t=0.00305; **0** nonincident boundary triangle intersection pairs in the geometric audit.
- Increment 36: converged at t=0.003055, in 54 Newton iterations, residual 0.0004248838594 < 0.0005; inner=0, outer=0; **0** boundary triangle intersection pairs.
- Increment 37: converged at t=0.003060, in 54 Newton iterations; inner=0, outer=0; **0** boundary triangle intersection pairs.
- Increment 38: entered Newton iteration 5 (no accepted state) before test execution limit; no active self-contact when step began.
- No separate physical second ring verified, and no production validation is claimed.

The tests above ran with the same `contact2.f` and projection implementation. After the last bounded run, an **outer-surface normal-sign selection correction** was added to `facet_search.f` and recompiled successfully, but the OUT-contact branch has not been numerically exercised.

## What actually changed

1. `contact2.f: SNODE_HT_IN` and `SNODE_HT_OUT`: call new `SELF_FACET_CANDIDATES` each increment rather than the original node-to-node candidate test, preserving separate inner/outer search sets.
2. `facet_search.f: SELF_FACET_CANDIDATES`: for each surface grid node, search nonlocal quadrilateral facets whose *expanded face bounding boxes* contain the node; verify their normal orientation, compute projected distance, and assign the closest facet's vertex as a seed for the legacy `SCSEGM_*` segment builder. Restricted to one axial-search direction to avoid duplicates in the legacy node numbering, currently excluding faces closer than seven axial rows in the structured grid. These are mesh-specific, experimentally chosen rules and can miss valid contacts.
3. `facet_search.f: SAFE_QUAD_PROJECT`: damped, clamped Gauss–Newton closest-point projection for four-node quadrilateral faces. Explicit step-tolerance and degeneracy protection replace the unbounded, one-sided `C_POINT` iteration.
4. `contact2.f: CONTACT_SEARCH`: scan available target *facets* and select by closest 3D projected point instead of choosing the nearest target vertex; use contact face coordinates and tangents at that projection for assembly. Faces coincident with a hitting node are skipped in the self-contact branch.
5. `facet_search.f`: a signed-gap candidate filter (`gap <= 5e-4`) stops the solver from constraining **open** noncontacting faces. This is a discovery filter, not a full unilateral active-set/complementarity update.
6. `INPUT_INC35`: conservative initial step reduced from 1.875e-5 to 5e-6. Force-residual tolerance remains unchanged (`5e-4`).
7. `build_sparse.sh` and `build_windows_ucrt64.sh`: compile/link the new `facet_search.f` source.

## Important unresolved implementation issues

- Original self-contact restart does not serialize the pair active set. Continue from pre-contact increment 35 *without stopping the executable* to retain state. Do not treat an increment-37 restart as equivalent to uninterrupted continuation.
- The `CHECK_CONTACT` routine still uses a multiplier/force-based heuristic rather than a rigorously checked gap and compressive reaction complementarity condition (`gap>=0`, multiplier>=0, gap*multiplier=0 with sign-appropriate convention).
- Legacy DOF renumbering assumes disjoint hitting/target node roles and may fail if a node participates in both; a general robust contact solver must use unique DOFs with independent constraint connectivity.
- The facet candidate finder currently chooses one target facet per hitting node and only searches nonlocal patches in one axial direction; it may fail at later folds and cannot prove self-contact is globally enforced.
- Accepted-step rejection/retry based on geometric intersections is not implemented **inside** the Fortran integrator. Use the independent `penetration_audit.py` before trusting stress contours.
- The geometric audit counts intersecting triangle pairs, not independent physical contacts; adjacent triangles sharing node IDs are excluded but further topology and numerical conditioning checks remain necessary.
- The high stress results inherited from checkpoint 35 remain unvalidated and must not be interpreted as material-response verification.

## Linux/WSL2 usage

Dependencies: `gfortran`, `gcc`, SuperLU development libraries (`-lsuperlu`), Python 3; optional `vtk`, NumPy and SciPy for intersection auditing.

```bash
bash build_sparse.sh
bash run_contact_trial.sh
```

The launcher backs up any `RESTART` folder, restores the saved precontact increment-35 checkpoint, and runs through later steps without an external execution limit. The original time limit in this test was outside the executable. Do not delete checkpoint folders unless making an intentional new start from zero.

## MSYS2 UCRT64 Windows usage

```bash
bash build_windows_ucrt64.sh
bash run_windows_ucrt64.sh
```

The Windows build script has been updated and checked for syntax, but **native Windows execution has not been tested here**. Keep the old working source separately.

## ParaView and geometry audit

```bash
python export_fem_outputs.py --root .
# On an environment with VTK/NumPy/SciPy installed:
python penetration_audit.py PARAVIEW_EXPORT/cylinder_inc0037.vtu
```

Open `PARAVIEW_EXPORT/cylinder_timesteps.pvd` in ParaView for the saved frames. The package contains actual accepted VTUs for 35, 36, 37. No failed-increment or second-ring geometry is presented as converged.

## Attribution

Original self-contact code written by Sarma, according to the project owner's attribution; the original file does not contain an explicit author header. This branch contains experimental modifications and is not proof of original author approval for redistribution.
