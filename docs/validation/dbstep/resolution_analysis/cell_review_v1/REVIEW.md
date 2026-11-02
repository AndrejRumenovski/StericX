# Independent review of the actual-cell diagnostic enclosure

Reviewed by the independent-method agent before the full convergence campaign.
The reviewed source hashes, command, seven passing focused tests and unchanged
source check are in `receipt.json` and `tests.log`.

The core derivation is sound for the declared finite, normal molecular-scale
arithmetic domain:

1. Exact rational midpoint cells partition a finite Cartesian box. Explicit
   endpoint checks guarantee that box contains the entire evaluation sphere.
   The three maximal node-to-cell half-spans give a valid Euclidean radius bound.
2. The subtract / square / three-term sum graph has the stated γ₅ distance bound.
   The two exact rational inequalities establish the sufficient radial envelope
   8u times radius. The campaign's nonzero-difference guard keeps squared terms
   normal, and the coordinate cap excludes overflow. Zero-radius ghosts remain
   empty by the declared contract.
3. A selected cell can disagree with continuous ball membership only within its
   cell-radius plus radial-roundoff shell. Union / intersection boundaries are
   contained in the union of the individual boundaries, so summed shell volumes
   are conservative even with overlapping atoms.
4. Every arbitrary change to a binary decision can alter weighted cell volume by
   at most the largest cell volume. Converting actual unequal cell weights into
   the observed constant-weight count uses the selected count times the maximum
   weight discrepancy; all width and weight calculations are rational.
5. Adding and subtracting Qₒ gives
   `|Vₛ Qₒ/Qₛ − Vₒ| ≤ Eₒ + (Qₒ/Qₛ) Eₛ` for Qₛ > 0. This needs neither a subset
   relation between numerator / denominator masks nor Eₛ < Qₛ. The separate
   reported-value difference accounts for normalization rounding and clamping.
6. Independent analytical lens intervals and adaptive box enclosures provide
   separate continuum evidence. These conservative cell bounds are diagnostic:
   broad bounds cannot establish fine accuracy, and charging observed mask
   defects cannot substitute for the separate implementation-validation lanes.

Three review clarifications were incorporated before the campaign: ordinary
floating error-interval arithmetic became rational with outward endpoints;
derived native percentage is explicitly labelled; and the reconstructed DBSTEP
fast denominator mask is identified as a count-equivalent auxiliary mask rather
than a directly observed per-node oracle. Equality to the upstream count suffices
for the count-estimator enclosure, while pointwise scientific claims require
separate direct evidence. No numerical residual fitted a constant or threshold.
