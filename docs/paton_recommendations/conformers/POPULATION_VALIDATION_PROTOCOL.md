# Independent population validation protocol

This protocol is written before running the independent population observer.
It tests arithmetic and metadata enforcement, not the physical validity of a
particular conformer search or electronic-structure method.

## Reference and acceptance rules

The reference is Python's standard-library Decimal arithmetic at 80 decimal
digits, evaluated directly as `g * exp(-(E - min(E)) * conversion / (R*T))`
followed by division by the sum. It does not import StericX or reproduce its
shifted-log implementation. Every numeric input is interpreted as the exact
binary64 number supplied to Rust (`Decimal.from_float`). Conversion factors are
derived independently from SI joules: exact Avogadro constant
6.02214076e23 mol^-1, exact Boltzmann constant 1.380649e-23 J/K, exact elementary
charge 1.602176634e-19 C, thermochemical kcal = 4184 J, and CODATA 2022 hartree
energy 4.359744722206e-18 J. Primary source:
[NIST CODATA table](https://physics.nist.gov/cuu/Constants/Table/allascii.txt).
The existing Rust R is retained; its binary64 rounding is tested against the
independent exact-SI reference, not copied into it.

For cases with at most 4096 states, each normalized weight and the sum must
agree within **32 binary64 epsilons absolute** (7.105427357601002e-15). This is
a numerical regression gate allowing rounding in conversion, subtraction,
division, logarithm, exponential and normalization; it is not a physical
accuracy claim or a rigorous universal libm bound. For nonzero weights at least
1e-300, additionally require **2e-12 relative error**: exponential sensitivity
scales with the magnitude of its exponent (up to about 700 in this range).
Neither criterion is fitted to observed residuals. Values below 1e-300 retain
the absolute gate; exact underflow behavior is tested separately. The Decimal
reference may replace an exponent below -10000 by zero, far below the smallest
binary64 subnormal, to avoid unbounded reference exponentiation.

All emitted weights must be finite and nonnegative. Tests retain input order,
raw requests and responses, reference weights, per-case errors, interpreter,
source and executable hashes, commands, and failures. Invalid input cases must
produce structured errors. A failing gate is reported and investigated; the
threshold is not loosened after seeing a failure.

## Cases fixed before observation

* Two states separated by RT ln(2), in each supported unit and at multiple T.
* Equal-energy multiplicities 1, 2, 4 and explicit enumerated duplicates.
* Multistate unequal energies, all four units, positive common offsets and
  negative electronic absolute energies; permutations preserve identity.
* All three energy kinds; Gibbs free energies at their evaluation temperature.
* Tiny positive temperatures, binary64 extreme energies/temperatures, large
  positive separations, exact underflow and a 4096-state ensemble.
* Missing energy metadata, mixed units, kinds and energy models, invalid
  temperatures, zero degeneracy, absent provenance and Gibbs T mismatch.

Closed-form Rust tests additionally check 2:1 two-state ratios, 1:2:4
degeneracy ratios, equal explicitly enumerated states and temperature behavior.
They are supplemental tests, not the independent reference. No energies in
this test are represented as historical experimental conformer energies.

The observer emits the exact input bit patterns accepted by Rust. The reference
uses those patterns, isolating population arithmetic from JSON parser rounding.
Any difference from the requested Python binary64 energy is separately recorded
and remains visible; it is not relabelled as a thermodynamic error or hidden by
rounding the reference. Invalid nonfinite JSON values are covered directly in
Rust because standard JSON cannot represent them.
