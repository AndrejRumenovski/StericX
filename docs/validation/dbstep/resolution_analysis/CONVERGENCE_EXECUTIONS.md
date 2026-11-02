# Convergence execution history

The accepted full campaign is
`.stericx/dbstep_resolution/convergence_v2`, with its standalone immutable
manifest, source snapshots, dependency/environment bindings and all raw
responses. It completed with exit code 0. The compact report and all 320
comparison rows are in
[convergence_v2/CONVERGENCE_REPORT.md](convergence_v2/CONVERGENCE_REPORT.md).

The first two-case smoke run completed 32 comparisons and exercised the raw
adapters and analytical/error-bound machinery. Its source is preserved verbatim
as `.stericx/dbstep_resolution/convergence_smoke_v1_controller.py`, matching the
recorded executing source hash. It was followed by full-grid native disputed-node
trace capture, outward error-interval endpoints and provenance improvements.
The smoke result does not substitute for the accepted full campaign.

The first full attempt (`convergence_v1`) stopped before scientific computation.
A controller guard incorrectly required the imported DBSTEP source path to equal
the checkout path. The pinned environment installs an identical source copy
under `site-packages`: both source SHA256 values were
`e3f1318eda340097aa21753e701498ac2bb497f8a973dcd12e298e27d9e7d48c`.
The corrected guard checks source-byte identity and binds both installed and
checkout paths. No scientific constant, numerical method or tolerance changed.
The failed log and exact controller remain preserved.

After the accepted campaign, the independent reference tests and convergence
derivation tests passed: **22 tests and 10 subtests**. Controllers pass Ruff
checks and formatting. The [execution receipt](convergence_execution_receipt.json)
binds the command, final logs, earlier attempts and summary/report receipt.

Reproduce the campaign in a new directory using the pinned reference environment:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 taskset -c 0,1 \
  .stericx/paton/reference/venv/bin/python \
  scripts/analyze_dbstep_volume_convergence.py \
  --observer .stericx/dbstep_resolution/mask_build_v2/bin/stericx-audit-observer \
  --output .stericx/dbstep_resolution/convergence_new_run
```

Regenerate the report from a completed, hash-verified campaign:

```sh
.venv/bin/python docs/validation/dbstep/resolution_analysis/summarize_convergence.py \
  --run .stericx/dbstep_resolution/convergence_new_run \
  --output docs/validation/dbstep/resolution_analysis/convergence_new_summary
```

These controllers report evidence; separate independent review decides admission
under the approved release policy. No release or publication claim follows solely
from completion of this campaign.
