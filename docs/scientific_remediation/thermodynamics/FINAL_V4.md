# Final v4 thermodynamics and kinetics replay

The frozen v4 executables and copied Python sources pass the scoped numerical
and invalid-input replay. The final [source/build binding](final_v4_build_binding.json)
connects the [replay manifest](final_v4_replay_manifest.json) to the complete
frozen build; [summary](final_v4_summary.json) records the measured results.
These are exact durable copies of the final receipts, with adjacent SHA-256
sidecars. Their artifact paths still refer to the original retained run at
`.stericx/scientific_remediation/thermodynamics/final_v4/`.

The replay covers **744 kinetic scalar comparisons**, **51 aggregate scalar or
rejection comparisons**, **14 CLI cases**, **6 CREST cases**, **11 MMFF cases**,
and **24 ee cases**. Normalization of all **1,566 Ni-hDA source rows** preserves
every published target and records the reaction condition as **353.15 K**.
The maximum normal-rate relative error is `5.2921798585206193e-8`; maximum ee
absolute error is `3.6928646096612283e-6` percentage points.

| Bound artifact | SHA-256 |
| --- | --- |
| Complete v4 build manifest | `d79fa7ec5c9e6c1254aa3bd1de665fe4338dde2aabe2401b7bf79a36b5a7bdd3` |
| Copied native executable | `102b6883c639bfc9ef211a88abae02391419cbd277c498098fa2bdfc7008206f` |
| Copied observation executable | `d6c5e20c8e012d1d9e5ff5afc5a8438697e0c1d11919b614ced69f0b8616ccc4` |
| Copied `prepare_data.py` | `6b9e49446c2b475a7093785ecb8993a29193dea0a0f1c26556275aeef7fe41f0` |
| Copied `stericx_quantum.py` | `df26cf4b9544d09d41b96e5827b9aff06ad2afbe562f4e000d6d235bbbb28c32` |
| Python observations | `d64c30ae8de95d41d71353debc4a7f823a383495a59d629cba65d3e64f038a06` |
| Scoped replay manifest | `e8260d0a7a5ceaa68dbd90a3283657c1a90e8cda0d6596526fafae8cffd36037` |
| Source/build binding | `dac4e878398726e782cab44cca1f803eae0702946307284398897bc66ec9b850` |

The copied executables match the complete v4 build exactly. All eight scoped
Rust/Python source files match its frozen source snapshot, and every replay
artifact was rehashed when publishing this summary. Python is 3.12.13.
Import-generated `.pyc` files are hashed runtime caches, not source evidence.

This result retains the limitations in [RESULTS.md](RESULTS.md) and the original
audit classifications. It does not establish conformer-sampling completeness,
free-energy accuracy, chemical calibration, or an absolute catalytic rate from
a barrier difference. Nine kinetic requests are invalid or exceed the bounded
Decimal comparison; focused checked-API tests separately cover extreme inputs.
External CREST execution was mocked for command/cache tests. MMFF status 1
remains an explicitly retained iteration-limit result, and missing populations
retain the documented uniform fallback. Correct reaction temperature metadata
does not imply historical conformer populations were recomputed at 353.15 K.
Unsigned ee cannot identify an absolute enantiomer, and source ligand 2064's
ee/ddG conflict remains unresolved. No historical report or evidence was
rewritten to create this receipt.
