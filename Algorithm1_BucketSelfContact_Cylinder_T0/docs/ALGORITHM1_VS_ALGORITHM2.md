# Contact-algorithm comparison

| Stage | Algorithm1 | Algorithm2 |
|---|---|---|
| Original method | Sarma's structured cylindrical `BUCKET` grid (`Z`–`theta`) and inner/outer node lists | Inherits solver; replaces candidate-finding stage with facet-aware searches |
| Search | `SNODE_HT_IN/OUT`, typically global every 20 increments or when no hits, otherwise local | `SELF_FACET_CANDIDATES` each increment, expanded facet bounding boxes |
| Target | Legacy node/segment proximity | Closest point to candidate quad face with bounded projection |
| Geometry screen | Legacy contact assembly + surface normals | Additional signed gap, normal-orientation and topological exclusion screen |
| Constraint mechanics | Legacy Lagrange-multiplier / active-set logic | **Still** retains that legacy enforcement; not a new penalty formulation |
| Validation | Early t=0 smoke test only, previous folding led to penetration | Early t=0 smoke test and a few small-step contact-prelude states only; full collapse NOT validated |

An experimental search-and-projection change is **not** proof the wall stops penetrating. Both branches require independent acceptance criteria. See respective README and known limitations.
