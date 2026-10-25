# Algorithm2 changes compared with earlier cylinder code

**Source provenance:** Sarma is identified by the project owner as the original self-contact author. An explicit Sarma author comment was not present in `contact2.f` in the original supplied archive. Preserve that attribution and seek permission before redistribution.

## Solver changes carried forward

1. `newrap1.f`: corrected skyline sentinel capacity, residual boundary-condition array limits, incomplete DOF masks, and added dynamic step growth/cutback and iteration control.
2. `contact2.f`: corrected earlier contact pair scratch array bounds and state persistence within a continuous run; experimental inner/outer contact node candidate logic.
3. `colsol_sparse.c` and `solvers.f` build: sparse SuperLU rather than duplicate legacy skyline COLSOL implementation.
4. `read.f`: initialise the `CYLINDER_RAIL` switch so mesh reading is deterministic.
5. `contact2.f`: experimental Algorithm2 facet-based candidate lookup and closest-contact-facet selection, with original multiplier constraint machinery still in place.
6. `facet_search.f`: expanded quad bounding-box tests, topological neighbour exclusion, orientation and signed-gap prefilter, clamped damped closest-point projection.
7. `scripts/`: fresh-start and optional pre-contact restart commands, no accidental checkpoint 35 import in the default run.
8. `visualization/`: ParaView VTU/PVD, fixed X-Z GIF, per-increment residual history and local no-cache dashboard, separate inner/outer ring candidate monitor and geometry intersection audit.

`FACET_CONTACT_CHANGES.patch` and `FACET_SEARCH_NEW_SOURCE.patch` are included under docs.
These patch files document the Algorithm2 facet trial against the preceding global-search branch, not the entire historical code lineage. Early source-change logs and exact original-to-current patch are available in earlier project packages; do not call these facet-only patches a complete original-to-current diff.

## Known issues

The source still uses a legacy multiplier active-set heuristic and mesh-specific search exclusions. It has no fully general crossed-facet contact enforcement, and no in-solver geometric rejection of accepted steps. It has not achieved independently verified nonpenetrating collapse or a second separate ring contact. The old checkpoint format does not serialize all active pair state, so mid-contact restart can change the physics. High stress values from old checkpoints remain suspect.
