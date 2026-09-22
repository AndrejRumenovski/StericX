# Planar clearance at input coordinate precision

Implemented after the full corrected-capture review and a retained failing-before regression.
This changes only the existing planar fallback (opposing unit-vector sum squared
at most 1e-4); that existing branch threshold is not fitted or changed here.

Each finite input f32 Cartesian component represents its rounding cell. Use half
of the larger adjacent f32 spacing as a conservative symmetric component
halfwidth (at the largest finite magnitude use the finite inward spacing).
For a donor component d and another atom component a, form the interval
[a-low - (d+high), a+high - (d-low)] in f64. For either existing one-angstrom
trial direction ±n, subtract that represented direction component, square each
component interval, then sum its three bounds. The lower bound is zero when the
component interval crosses zero. Round arithmetic bounds outwards to adjacent
f64 values. The minimum of the atom-wise lower bounds and minimum of their upper
bounds encloses the nearest-heavy-atom clearance for that trial direction.

Select +n only if its lower clearance bound exceeds the -n upper bound; select
-n only for the converse. Overlap means the sign is unresolved at input
coordinate precision and requires an explicit center. These intervals do not
claim experimental coordinate uncertainties or a physical lone pair. They state
when the existing geometric trial-direction model may resolve its sign from
rounded input coordinates. The trial normal is held at its computed represented
value for the comparison. This is not a tolerance applied to descriptor output.

Read-only diagnostic: all 124 original planar base/rotation/permutation requests
have overlapping coordinate-derived intervals. The previous exact-equality test
rejected 29 but selected a sign for 95. The largest diagnostic squared-distance
interval width was 5.38839e-5 Å², a measured bound arising from those input floats,
not a chosen acceptance threshold. A nearby off-plane carbon at 2.5 Å provides a
separate unambiguously resolved control.

The implemented comparison uses outward-rounded f64 primitive operations. All
124 original variants now return the explicit ambiguity error; a clearly
off-plane obstruction retains the correct resolved sign. Exact conditional
interval endpoints and the independent arithmetic program are retained in
`outward_interval_evidence.json` and `diagnose_intervals.py`. Inferred-normal
uncertainty is not propagated: these are conditional coordinate-rounding bounds,
not a rigorous enclosure of every possible inferred geometric frame.
