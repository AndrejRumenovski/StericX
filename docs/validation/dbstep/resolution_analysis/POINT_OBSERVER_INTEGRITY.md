# Native point-decision observer integrity

The new [observation adapter](../resolution_observer/RESOLUTION_OBSERVER.md)
captures actual native point decisions rather than reconstructing the native
distance predicate in Python. It is the system under investigation, not an
independent reference. Production code, existing observers, historical evidence,
and the frozen baseline are unchanged.

The accepted build uses `.stericx/dbstep_resolution/mask_build_v2`, with its own
source copy and Cargo target. The complete production source set matches the
frozen pre-investigation baseline. The compiled observation source begins with
the unchanged 65,574-byte production volume source; its SHA256 is
`6e4e0f368d9acf8f4b20327d5281384a42836d56e039ce4824571a2c45ec54bd`.
The observer executable SHA256 is
`40c7953b3ec9253247613a5f76b436601b99a567a58b65b642a3a228965f00a0`.
The [exact adapter diff](mask_observer_adapter.diff) shows the old source files
are unchanged apart from the explicit include and operation dispatch additions.
Added files include the observation append, documentation and build/check
controllers. [Build provenance](mask_observer_build_receipt.json) binds the diff,
source hashes, toolchain-backed build manifest, command and executable.

For original failure `SIG-NIHDA-401/0/sphere_radius=2.5`, the new observer retains
5,565 generated sphere nodes in all three native coordinate frames:

| Native plane atom | Occupied nodes | Aggregate/regional reconstruction |
|---|---:|---|
| 15 | 1,396 | Bit-for-bit equal |
| 13 | 1,398 | Bit-for-bit equal |
| 10 | 1,397 | Bit-for-bit equal |

Each occupancy mask entry calls the unchanged private `occupied_volumes` with a
one-point slice. The full native aggregate is called independently on the entire
point set; equality is therefore checked rather than assumed. Counting masks
and native sign-based region IDs reconstructs total, quadrant, octant, near and
far volumes with the exact f32 native normalization. Native generated-grid and
supplied-point operations yield identical masks and outputs. Legacy geometry
and raw-occupancy operations match the frozen baseline observer and original
historical frame values. Empty point sets, invalid density and underflowing
sphere volume produce explicit errors and no masks.

[Integrity results](mask_observer_integrity_result.json) and the
[evidence receipt](mask_observer_evidence_receipt.json) bind all raw requests,
responses and errors. The check is reproducible with:

```sh
.venv/bin/python docs/validation/dbstep/resolution_observer/check_observer.py \
  --baseline-binary .stericx/dbstep_resolution/baseline/build/bin/stericx-audit-observer \
  --binary .stericx/dbstep_resolution/mask_build_v2/bin/stericx-audit-observer \
  --output .stericx/dbstep_resolution/mask_integrity_new_run
```

The first build attempt at `mask_build` compiled both binaries but its immutable
input guard rejected an unrelated replay-controller change during compilation.
That attempt, its error and diagnostic integrity run remain preserved. It is
excluded from accepted evidence. The accepted v2 build was run after the replay
controller author held edits and all guards passed. This is an evidence-integrity
failure, not a scientific result or an implementation correction.

These checks establish trustworthy observation of the native kernel. Independent
finite-grid and continuum comparisons remain separate scientific validation
requirements; this document does not claim those gates have passed.
