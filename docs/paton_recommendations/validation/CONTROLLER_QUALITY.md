# Validation controller quality checks

Maintained Paton study controllers now pass repository formatting and lint rules.
Their numerical equations, scientific constants, cohorts, thresholds and result
criteria are unchanged. Corresponding input/output arrays use `zip(strict=True)`
to reject length mismatches instead of truncating a comparison. Local variable
renames, import cleanup, formatting and executable modes complete the cleanup.
Exact provenance strings remain unchanged; a few long literal lines carry scoped
line-length annotations.

Frozen failed-run controllers in `docs/validation/dbstep/controllers/` retain
their recorded bytes and are excluded through a local Ruff configuration, as
existing sealed audits already do. Original executed controller versions remain
available in their study commits; evidence receipts continue to identify those
versions. This cleanup does not retroactively alter any archived result.

Ruff check and format verification pass for the maintained repository (233 files
at this checkpoint), excluding user-owned untracked `docs/media/`. Those user
files were inspected by a broad preflight but were not modified. Four focused
controller tests pass. Complete final scientific replay and engineering tests
remain separate gates.
