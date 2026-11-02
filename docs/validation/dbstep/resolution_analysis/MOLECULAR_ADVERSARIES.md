# Molecular transformations and focused regressions

The seeded campaign contains **49 requests**: PH3 and primary, secondary and
tertiary synthetic phosphines with and without hydrogen occupancy, plus the
original ligand-401 R=2.5 witness. Each has identity, randomized atom permutation,
signed-axis rotation, general rigid rotation and translation variants. Four
nearly collinear constructions complete the matrix.

**48 requests succeed; the exactly collinear request returns a structured
geometry error.** All accepted cases have identical individual occupancy masks
in Rust, DBSTEP and the independent binary64 reference. There are no unexplained
point decisions or count-to-volume arithmetic errors. Atom permutations preserve
all public volume-summary bits. The tested rotations/translations change no
occupied nodes; aligned-coordinate displacement is nevertheless retained for
each frame rather than claiming universal bitwise rotational invariance.

The first controller incorrectly expected a nonzero transverse offset of
`2^-24 Å` to be rejected. That expectation was wrong: with the specified axis,
the exact plane-normal magnitude is `2.25*abs(offset)`, which is nonzero. The
existing audited Rust fallback uses a binary64 cross product to resolve this
valid small plane. The corrected expectation rejects only zero offset. Both
attempts, original controller bytes, inputs, results and the explicit correction
remain in [the evidence directory](molecular_adversaries_v2/). No production
threshold or implementation changed to make this test pass.

Six focused integration regressions additionally verify:

* the original 32-node sphere reclip and six occupied-node changes per frame;
* exact and adjacent-float atom boundaries, with rational roundoff witnesses;
* rejection by the checker of a deliberately wrong stable kernel decision;
* quadrant/octant assignments on zero-coordinate planes;
* row-cache behavior under shuffled points and atom order;
* an empty atomic union preserving the sample domain with zero occupation.

All six pass using the unmodified frozen DBSTEP environment and the source-bound
Rust observer. Environments missing either prerequisite explicitly skip these
integration tests; such a skip is not accepted as reference validation.

The additional [public-summary review](point_domains_review_v2/public_summaries.json)
checks **138** mean-volume and percentage results from the original 69 cases.
Every result satisfies its unchanged original arithmetic tolerance after using
the common sample domain. The earlier frame-level evidence remains unchanged.

```bash
taskset -c 0,1 .stericx/paton/reference/venv/bin/python \
  scripts/analyze_dbstep_molecular_adversaries.py \
  --observer .stericx/dbstep_resolution/mask_build_v2/bin/stericx-audit-observer \
  --output .stericx/dbstep_resolution/molecular_adversaries_fresh
taskset -c 0,1 .stericx/paton/reference/venv/bin/python -m unittest discover \
  -s tests -p test_dbstep_point_domains.py -v
```

These results address the tested geometries. They do not turn f32 input geometry
into exact real coordinates or remove the resolution dependence of lattice
integration. Continuous-volume bounds and the final approved-policy review are
separate evidence.
