# Releasing StericX and minting a Zenodo DOI

A reusable runbook: substitute the version being cut for `X.Y.Z` throughout.
v0.1.0 and v0.2.0 were both cut this way.

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

- `Cargo.toml` — `version = "X.Y.Z"` (and `cargo build` so `Cargo.lock` follows).
- `CITATION.cff` — `version: X.Y.Z`, `date-released` = today.
- `CHANGELOG.md` — the top section is `## [X.Y.Z] - <date>` (not `[Unreleased]`),
  with one `### Added` / `### Changed` per section rather than duplicates.
- `.zenodo.json` — `version` bumped and the description covering what is new.
- `.zenodo.json` — title, description, creators, and `related_identifiers` current.

## 1. **(you)** Connect the repository to Zenodo — first release only

1. Sign in at <https://zenodo.org> with GitHub (Sign in → GitHub), authorizing
   the Zenodo GitHub app.
2. Go to <https://zenodo.org/account/settings/github/> and flip the toggle for
   `AndrejRumenovski/StericX` **On**. (If it is not listed, use *Sync* / re-check
   permissions.) Zenodo now watches the repo for new releases and reads
   `.zenodo.json` for metadata. **Already enabled** — subsequent releases skip
   this step entirely.

## 2. Tag and publish the GitHub release

The tag is what Zenodo archives. Create it from the release commit on `main`:

```bash
# annotated tag
git tag -a vX.Y.Z -m "StericX vX.Y.Z"
git push origin vX.Y.Z

# publish the GitHub release (notes from the matching CHANGELOG section)
sed -n '/## \[X.Y.Z\]/,/^## \[/p' CHANGELOG.md | sed '$d' | tail -n +2 > /tmp/notes.md
gh release create vX.Y.Z \
  --title "StericX vX.Y.Z" \
  --notes-file /tmp/notes.md
```

Write the notes to a file rather than using a process substitution: `gh` needs a
real seekable file, and `tail -n +2` drops the duplicated heading.

Publishing the release fires the Zenodo webhook; Zenodo mints the DOI within a
minute or two.

## 3. Record the DOI

1. Zenodo mints within a minute or two. The version DOI can be read straight from
   the public API without a browser:

   ```bash
   curl -s "https://zenodo.org/api/records?q=conceptrecid:21726666&all_versions=true&size=10" \
     | python3 -c "import sys,json;[print(h['metadata'].get('version'), h['metadata'].get('doi')) for h in json.load(sys.stdin)['hits']['hits']]"
   ```

   Zenodo issues **two** DOIs: a *concept* DOI (all versions, constant —
   `10.5281/zenodo.21726666`) and a *version* DOI for this exact release.
2. Put the concept-DOI badge at the top of `README.md` and replace the
   `ARCHIVE-DOI` placeholder in the "How to cite" section with the concept DOI:

   ```markdown
   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
   ```

3. Record the DOIs. Keep `CITATION.cff`'s top-level `doi:` as the **concept** DOI —
   a version DOI there is stale the moment the next release is archived, and the
   snapshot Zenodo takes at the tag would preserve the stale value. List every
   version DOI under `identifiers:` instead, and note the new one in the
   CHANGELOG section for that release. Commit as `docs: record Zenodo DOI`.

## Notes

- Re-releasing (e.g. v0.1.1) repeats steps 2–3; the concept DOI stays constant and
  a new version DOI is minted automatically.
- The archived artifact is a snapshot of the repository at the tag. Large
  gitignored caches (`.stericx/`, `data/external/`) are **not** included, by
  design — they are third-party or regenerable, and the copyrighted Science SI is
  never redistributed.
