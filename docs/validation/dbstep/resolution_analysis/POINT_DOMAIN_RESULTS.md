# Common-point-domain investigation

The corrected comparison supplies each captured native sphere grid directly to
unmodified DBSTEP atomic occupancy and an independently written union predicate.
It does not clip the same nodes again before comparing atomic decisions.
This changes a research comparison, not the Rust descriptor engine.

All **69 original cases / 207 orientations** have exactly equal individual
occupancy masks in Rust, DBSTEP and the independent direct reference. Native
integer counts reconstruct all aggregate, quadrant and octant output bits.
All common-domain absolute volumes satisfy the **unchanged original arithmetic
tolerance**. There are zero unexplained atomic or sphere decisions.

The historical reclip lane still fails in the same three ligand-401 R=2.5 frames:

| Plane | Rust occupied / total | DBSTEP reclip occupied / total | Original residual, Å³ |
|---|---|---|---:|
| 15 | 1396 / 5565 | 1390 / 5533 | -0.0239801436431 |
| 13 | 1398 / 5565 | 1392 / 5533 | -0.0241167166704 |
| 10 | 1397 / 5565 | 1391 / 5533 | -0.0240484301567 |

The same 32 nodes are removed from each denominator; six occupied nodes are
removed per frame. Exact rational margins and explicit binary32 operation traces
explain every affected sphere point. Their mathematical positions lie just
outside the exact sphere, while binary32 squared-distance arithmetic rounds
onto its inclusive boundary. The mathematical continuum boundary has zero
measure; these sampled nodes have finite estimator weight.

This is a **reference comparison error** in the asserted common-domain equality,
with an underlying **finite-precision integration difference**. No atomic
occupancy defect in StericX or DBSTEP was found in these 207 frames. The separate
DBSTEP denominator/numerator inconsistency documented in
[DBSTEP_SOURCE_REVIEW.md](DBSTEP_SOURCE_REVIEW.md) is not the cause of this
original tree-path failure.

The [complete frame records](point_domains_v1/frames.json) retain every case,
accepted frame geometry/radius, residual and rational boundary witness.
The [receipt](point_domains_v1/receipt.json) rehashes all 213 raw artifacts
(approximately 313 MiB); raw requests and responses remain in
`.stericx/dbstep_resolution/point_domains_v1`. Existing input reports are
independently hash-checked by the unchanged case loader.

Reproduction in the frozen DBSTEP environment:

```bash
taskset -c 0,1 .stericx/paton/reference/venv/bin/python \
  scripts/analyze_dbstep_point_domains.py \
  --observer .stericx/dbstep_resolution/mask_build_v2/bin/stericx-audit-observer \
  --output .stericx/dbstep_resolution/point_domains_fresh
.stericx/paton/reference/venv/bin/python scripts/report_dbstep_point_domains.py \
  --input .stericx/dbstep_resolution/point_domains_fresh \
  --output /tmp/stericx-point-domain-report
```

These are investigation results, **not a revised release-gate pass**. Analytical
accuracy, adversarial geometry, convergence, complete replay, packaging and
matched performance remain separate requirements.
