# Changes in adaptive cylinder continuation — 24 October 2026

## Inherited source

All prior code fixes are described in `CHANGELOG.md`. The last independently strict-tolerance checkpoints were increments 31–33; these were produced after the original t=0 calculation converged to increment 30. The revised code is for development/diagnostic use and not independently validated against a reference FEM impact result.

## New edits made in this continuation

1. `newrap1.f`, `STRESS_ANALYSES`: Retain step size between successive iterations in the **same process** rather than resetting `DT=FDT` each time; store and report planned `FDT`.
2. `newrap1.f`: After two accepted steps using at most 60 Newton iterations, increase planned next `DT` by ×1.25 (upper cap 1e-4); after a slow accepted step (≥80 iterations), reduce ×0.8; very slow (≥100), halve it. Maximum Newton iterations in input = 200. The convergence tolerance is untouched at 5e-4.
3. `newrap1.f`, `NEW_RAP`: At maximum Newton iterations without convergence, attempt a half-size time-step cut, retaining current uncommitted state and returning to the existing trial-initialisation label. The limit is seven cuts, minimum 1e-7. **This retry path was added but not yet validated against a known reference case**; it may need more complete rollback of auxiliary contact history.
4. `newrap1.f`: If either inner or outer self-contact starts, reduce the *current* time step to 1e-5 before factorisation if it was larger, and cap next step at 1e-5 whenever self-contact is active. Self-contact detectors remain independent.
5. `contact2.f`, `SCSEGM_IN` and `SCSEGM_OUT`: The list of **unique target vertices** was incorrectly declared as `4*SIZE_HIT`, whereas up to four times as many distinct vertices as facets may be required. The actual caller allocates the full 61×24 node-grid capacity. The dummy array dimensions now match it, with explicit overflow STOP and comparison against only inserted vertices (`IDX2-1`). This removes the overrun first seen in increment 39.
6. `contact2.f`, `SELF_CON_STAT`: `DCNDN_IN` and `DCNDN_OUT` dummy extents previously depended on stale `SIZE1+SIZE2` and `SIZE3+SIZE4` values that are updated inside the routine. They are now dimensioned as `2*SELFCON_NODES_IN/OUT`. `newrap1.f` caller workspaces are enlarged correspondingly. This fixes the runtime bounds overrun first observed at increment 40.
7. `adaptive_dashboard.py`: Generates a refreshed local HTML dashboard, a labelled Newton residual plot, accepted/rejected iterations, contact activation events, step increases/reductions, retry markers, and potential axial rings. Ring detection based on contact-search rows **does not independently verify a second physical ring**.
8. `run_with_monitor.sh`: Backs up the saved checkpoint, uses the current Linux solver executable, monitors its stdout, and attempts to restore the planned next time-step from matching log metadata on a new program launch. It does **not** restore the complete inner/outer active pair lists; avoid unnecessary process restarts while self-contact is active.

## Evidence, run history and reservations

See `RESULTS_ADAPTIVE.md`, the log files, CSV summaries, and `monitor.html`. A solver diagnostic that identifies contact does not imply a converged accepted increment. Only increments for which `NEWRAP CONVERGES` and strict residual < 5e-4 agree count as accepted.

The attempt following increment 38 found the **first inner self-contact set** (four hitting nodes at axial row 18) in increment 39, and a follow-up attempt strictly converged that step. At increment 40, nine inner hitting nodes were identified, but the run exposed the stale `DCNDN_IN` declaration described above. These source corrections compile with `-fcheck=bounds`. The continuous-from-38 retest of both corrections should be consulted for any later acceptance claims.

Despite two tolerated restart continuations to reach this stage, the entire history is **not** a single uninterrupted t=0 run. The original restart active-contact history is incomplete. Before engineering sign-off, verify signed normal gaps, force direction, contact pair legality and energy/momentum balance through the full crushing history. Neither outer contact nor a second physical contact ring has been independently verified yet.

## Additional corrections after increment-40 diagnostic

9. **Persistent target-node lists:** `S_TNODE_IN` and `S_TNODE_OUT` were local, temporary arrays in `SELF_CON_STAT`, while `S_HNODE_IN` and `S_HNODE_OUT` were held by `STRESS_ANALYSES`. This asymmetry caused undefined previous target IDs on the next increment. Both target lists are now declared in `STRESS_ANALYSES`, initialised once, and passed explicitly into `SELF_CON_STAT` and onward into the local searches. These arrays **persist across successful increments within a single executable run**, but remain absent from the original restart file format. Thus mid-contact restart is still not equivalent to an uninterrupted calculation.
10. **`EXHIT` storage:** `EXHIT_IN/EYHIT_IN/EZHIT_IN` and the outer equivalents in `SELF_CON_STAT` were declared using the previous increment's `SIZE1/SIZE3`, even though the hit count may increase when current contact is detected. Their dummy dimensions now reflect the full original inner/outer contact surface-node capacities, consistent with the parent declarations.
11. Recompiled with `gfortran -O2 -fcheck=bounds -fbacktrace` and SuperLU. Source changes validated for compilation; downstream contact-tangent and physical gap/force consistency remain unvalidated.
