# CI maintenance after campaign completion

The two pushes ending at `6933b5bc` and its reworded equivalent `c4aefa6a`
failed the Python lint step with 1,067 findings in the completed campaign's
research scripts. Rust formatting, Clippy, documentation and tests passed;
macOS/Windows Rust jobs and the Python wheel/source-install matrix passed.

The executed scripts are scientific evidence: their exact bytes are bound into
admission, benchmark and completion receipts. Formatting those archived scripts
would invalidate their recorded hashes. This maintenance change follows the
existing `docs/ultimate_engine` archive policy:

- `published_archive.json` records all 95 campaign files present at `c4aefa6a`,
  including the completion audit and all 36 Python scripts. Its SHA-256 is pinned
  in `scripts/verify_campaign_archive.py`.
- CI checks sizes, hashes, Python syntax and Git tracking. Missing files,
  changed seals, tampering and symlink replacements fail.
- `.ruff.toml` excludes only the 36 exact sealed Python filenames. CI requires
  that list to match the archive and rejects exclusions that would hide a newly
  tracked nested helper. New helpers inherit the repository's full Ruff rules.
- Both campaigns' archive guards have independent tamper, missing-file,
  path-safety, syntax and exclusion tests. Production Rust and Python binding
  checks remain unchanged.

This verifies published archive integrity. It does not rerun the full research
campaign or claim to verify workstation-only raw measurements in a bare clone.
No historical file, scientific source, acceptance tolerance or result was edited.
The root README now reports the latest controlled baseline/results and the
explicit certified API with its assumptions and unmet tolerance targets.

Run the maintenance guards from the repository root:

```sh
uv run python scripts/verify_campaign_archive.py
uv run python -m unittest scripts/test_verify_campaign_archive.py
uv run ruff check .
uv run ruff format --check .
```
