# Releasing StericX (v0.1.0) and minting a Zenodo DOI

This is the runbook for cutting a tagged GitHub release and archiving it to Zenodo
for a citable, immutable DOI. Steps marked **(you)** require account access and a
browser and cannot be automated from the repo; the rest are shell commands.

**Order matters.** Zenodo only archives a GitHub release if the repository's
Zenodo integration is switched on *before* the release is published. Do the steps
in order.

## 0. Pre-flight (all green already, but re-check on the release commit)

```bash
cargo test --all-targets
cargo fmt --check && cargo clippy --all-targets -- -D warnings
uv run ruff check . && uv run ruff format --check .
uv run python -m unittest tests/test_quantum_backend.py
```

Confirm the metadata is consistent:

- `CITATION.cff` — `version: 0.1.0`, `date-released` = today.
- `CHANGELOG.md` — the top section is `## [0.1.0] - <date>` (not `[Unreleased]`).
- `.zenodo.json` — title, description, creators, and `related_identifiers` current.

## 1. **(you)** Connect the repository to Zenodo — do this first

1. Sign in at <https://zenodo.org> with GitHub (Sign in → GitHub), authorizing
   the Zenodo GitHub app.
2. Go to <https://zenodo.org/account/settings/github/> and flip the toggle for
   `AndrejRumenovski/StericX` **On**. (If it is not listed, use *Sync* / re-check
   permissions.) Zenodo now watches the repo for new releases and reads
   `.zenodo.json` for metadata.

## 2. Tag and publish the GitHub release

The tag is what Zenodo archives. Create it from the release commit on `main`:

```bash
# annotated tag
git tag -a v0.1.0 -m "StericX v0.1.0"
git push origin v0.1.0

# publish the GitHub release (notes from the CHANGELOG 0.1.0 section)
gh release create v0.1.0 \
  --title "StericX v0.1.0" \
  --notes-file <(sed -n '/## \[0.1.0\]/,/## \[/p' CHANGELOG.md | sed '$d')
```

Publishing the release fires the Zenodo webhook; Zenodo mints the DOI within a
minute or two.

## 3. Record the DOI

1. **(you)** Open <https://zenodo.org/account/settings/github/> → the StericX
   entry now shows a DOI badge. Zenodo issues **two** DOIs: a *concept* DOI (cite
   all versions) and a *version* DOI (this exact release). Copy both.
2. Put the concept-DOI badge at the top of `README.md` and replace the
   `ARCHIVE-DOI` placeholder in the "How to cite" section with the concept DOI:

   ```markdown
   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
   ```

3. Add the DOI to `CITATION.cff` (top-level `doi:` = the version DOI) and, if
   desired, to `.zenodo.json`. Commit as a small `docs: add Zenodo DOI` change.
   The `.zenodo.json` used for the *next* release will then carry it forward.

## Notes

- Re-releasing (e.g. v0.1.1) repeats steps 2–3; the concept DOI stays constant and
  a new version DOI is minted automatically.
- The archived artifact is a snapshot of the repository at the tag. Large
  gitignored caches (`.stericx/`, `data/external/`) are **not** included, by
  design — they are third-party or regenerable, and the copyrighted Science SI is
  never redistributed.
