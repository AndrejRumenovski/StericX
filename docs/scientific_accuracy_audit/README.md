# Scientific accuracy audit

Read the [main report](SCIENTIFIC_ACCURACY_AUDIT.md), the [individual claims and classifications](CLAIMS.md), and the [reproduction instructions](REPRODUCE.md). The audit records observations and independent comparisons against StericX commit `6393aafe0d983e504baf8abc1e18dc2a0f0d40e7`; publishing these artifacts does not change the scientific implementation or resolve the findings.

## Restore and verify a GitHub clone

From the repository root, using Python 3:

```sh
python3 docs/scientific_accuracy_audit/publication/manage.py restore
python3 docs/scientific_accuracy_audit/publication/manage.py verify
```

Reports, scripts and selected metadata can be read directly in GitHub. The larger evidence collection is stored as ordered parts of one losslessly compressed tar archive, each at most 48 MiB. Restoration checks every part and the combined archive against `publication/manifest.json`, then restores the listed files to their original audit-relative paths and verifies their sizes and SHA-256 hashes. Inputs, native outputs, independent reference results and recorded failures retain their original bytes. Frozen executables retain their executable modes. Restoration is repeatable and refuses to replace an existing file with different bytes.

The new verifier checks the published files, the original sealed manifest, and the current repository source/input files recorded in `manifest_initial.json`. An audit-only descendant commit is valid: its Git commit ID need not equal the frozen commit if those recorded repository files remain identical. The verifier reports its scope as **published evidence** and exits unsuccessfully for a missing or changed published file or a changed recorded repository file.

Journal articles and other third-party source captures excluded under the repository's source-publication policy remain local. Their acquisition metadata, original paths, sizes and hashes are retained; `publication/manifest.json` lists each exclusion and its reason under `local_only`. The verifier explicitly lists absent local-only captures and verifies their hashes when they are present. Passing published-evidence verification does **not** mean the entire original source corpus is included in this clone. Reacquisition may require access to the original provider, and a newer response may not have the frozen hash.

For a scratch extraction, copy the published audit directory, including its manifest and bundle parts, to the scratch location, then use:

```sh
python3 docs/scientific_accuracy_audit/publication/manage.py restore --audit-dir /path/to/scratch/audit
python3 docs/scientific_accuracy_audit/publication/manage.py verify --audit-dir /path/to/scratch/audit --repo-dir "$PWD"
```

## Original full verification and replay

The sealed [REPRODUCE.md](REPRODUCE.md), `scripts/audit_package.py`, and `rechecks/recheck_current_system.py` retain their original assumptions and bytes. Those original verifiers require the original Git commit, so use the new publication helper above on a current GitHub clone.

Historical audit snapshots are excluded from style rewriting to preserve their hashes; publication helpers and production code remain linted.

To perform the original full verification, create a separate worktree at `6393aafe0d983e504baf8abc1e18dc2a0f0d40e7`, copy the restored audit into it, and supply the locally retained or reacquired source captures at their recorded paths with matching hashes. Then follow [REPRODUCE.md](REPRODUCE.md) there. The original scripts are deliberately unchanged; a publication verification pass does not stand in for their complete checks. Replaying observations checks reproducibility, while recomputing independent references requires the recorded tools and inputs. Keep new outputs in a scratch worktree or a new experiment directory so the frozen evidence remains unchanged.
