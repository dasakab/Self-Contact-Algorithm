# Cylinder_250 heap-corruption diagnostics — 24 Oct 2026

## Base
Copied source from `cylinder_250_adaptive_continuation` (latest contact/target-list and adaptive-time-step modifications), and reproduced from the independent pre-contact increment 38 restart (`t=0.00306125`). This is not a fresh uninterrupted t=0 run.

## Confirmed defects and minimal source edits

### H1 — `newrap1.f` skyline LOC sentinel overflow
AddressSanitizer reported a 4-byte heap-buffer-overflow WRITE in `ADDRESS`, `skyline.f:53`, through `STRESS_ANALYSES/newrap1.f:649` before increment 39's first nonlinear solve. `ADDRESS` writes `LOC(ND+1)`, but the original `LOCA(NTDOF_A)` / `LOCB(NTDOF_B)` could be one integer short when all DOFs are active.

Fix in the parent `STRESS_ANALYSES`: `INTEGER LOCA(NTDOF_A+1),JCHT_A(NTDOF_A)` and similarly for `LOCB` with `NTDOF_B+1`. The ADDRESS routine and contact formulation are unchanged. Exact diff: `HEAP_ADDRESS_FIX.patch`.

### H2 — `UNBALANCE_FORCE` reads beyond contact-condition buffers
Once H1 was corrected, AddressSanitizer reported a 4-byte heap-buffer-overflow READ in `UNBALANCE_FORCE`, `newrap1.f:2895`. It looped `K1=1,IA` but read `CDISP_A(K1)`, allocated only for `3*NCOUNT_A` ordinary contact DOFs (similarly for body B). With self contact, IA includes other DOFs.

Fix: add `NCOUNT_A,NCOUNT_B` to the routine arguments and make CDISP dummy sizes `3*NCOUNT_A` and `3*NCOUNT_B`. Limit masking loops to these actual extents; do not blindly treat self-contact components as CDISP flags. Exact diff: `RESIDUAL_MASK_FIX.patch`. This may change the residual/masking compared with prior runs, so the contact forces and remaining BC logic need independent verification.

### H3 — `SELF_CON_STAT` segment workspace sized by stale segment count
A contact-routines-only AddressSanitizer build reproduced a heap-buffer-overflow WRITE in `SCSEGM_IN`, `contact2.f:3221` (`S_CONSEG_IN = 0`), called from `SELF_CON_STAT`. The caller declared `S_CONSEG_IN(4,SIZE_SEG_IN)` based on the **previous increment**. The callee was filling `4*SIZE_HIT_IN` facets using the new inner contact set, overflowing the stale allocation as hits increased. The same issue was present for outer contact, and for additional locally declared segment arrays.

Fix both inner and outer `S_CONSEG` workspaces to `(4,4*SELFCON_NODES_IN/OUT)` and `NCSEG`/`CONSEG` to the corresponding full surface capacities, while using `SIZE_SEG_*` only as the updated *active* count. Exact diff: `HEAP_SEGMENT_CAPACITY_FIX.patch`. This fix has been compiled; post-fix numerical continuation is tracked separately.

## Test evidence
- AddressSanitizer built with `-g -fcheck=all -fsanitize=address` (see `build_asan.sh` in the companion diagnostic folder). It first reported H1, and reported H2 on the immediately following test after H1 was fixed.
- Optimised bounds-checked SuperLU run with H1+H2 corrected: increment 39 converged with 52 Newton iterations, residual 0.000480624934789, tolerance 0.000500000023749, time 0.00307125. At increment 39 the inner detector saw 4 hitting nodes; outer 0.
- Increment 40 identified 10 inner hitting nodes, 28 inner contact segments, and 0 outer hits but still aborted during `CON_SEGMENT` allocation (`malloc(): unaligned tcache chunk detected`) in the H1+H2-only executable. The subsequent targeted sanitizer run identified H3. This log is retained as a before-H3 baseline; see `run_fixed_seg.log` for the after-H3 run.
- An instrumented restart *from increment 39* followed a different active-set path due to incomplete self-contact restart history, so it is for diagnostics only.

## State of numerical validation — updated after H3 correction

The H1/H2-only run still aborted at increment 40. A contact-routines AddressSanitizer trace confirmed H3 at `SCSEGM_IN` (line 3221), and the release build with H1+H2+H3 cleared that heap abort and executed the first increment-40 factorisation. With `DT=1e-5`, it did not reach force convergence before the bounded run was stopped. A controlled, **continuous from precontact increment 38**, smaller-step `DT=5e-6` test converged increments 39, 40, 41, 42, 43 and 44, each at the unchanged residual tolerance of 5e-4 (see `README_HEAP_CORRECTIONS.md` for exact residuals). Inner contact was active; outer contact and any separate second ring were not detected. The residual-mask change from H2 means comparisons with the old uncorrected run are not mathematically equivalent.

Only memory bugs H1, H2 and H3 are confirmed and changed by this study. Other legacy self-contact code may contain further errors. Element energy, contact-force equilibrium, gap signs, contact normal orientation, and independent time-step convergence are unvalidated. Mid-contact checkpoint directories are intended for viewing only since the archived restart format does not preserve the complete in-memory contact pair state.
