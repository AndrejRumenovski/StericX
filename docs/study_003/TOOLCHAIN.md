# Study 003 Quantum Toolchain

## Reproducible executables

`install_quantum_tools.sh` downloads the upstream binary releases, verifies
their archives before extraction, and installs them below the ignored
`.stericx/tools/` directory.

| Tool | Release | Release archive SHA-256 | Installed executable SHA-256 |
|---|---:|---|---|
| xTB | 6.4.0 | `c31f1c446a5a78a1e5e558b6e688904ae9b0398272b07f260f6e68a18fa27412` | `75e52241ed8a59fd1ac72f471e483ecb430754d52977a57fb855a0fb8a6d0863` |
| CREST | 2.12 | `c55e0f075a6223317b33a5f0fae593ce0ad55c1229c382937b0a0c2dcaf72ef6` | `ac9cbe9ea270aa2288d5d83b3c40400ad9c39124459b57d257e48c02217a5a6d` |

These versions reproduce the versions named by the public Kraken workflow.
The source conventions are the
[official Kraken repository](https://github.com/SigmanGroup/kraken),
[CREST documentation](https://crest-lab.github.io/crest-docs/), and
[xTB command-line documentation](https://xtb-docs.readthedocs.io/en/latest/commandline.html).

## Production calculation profile

The CREST ensemble command is equivalent to:

```text
crest input.xyz --gbsa toluene -metac -nozs -T THREADS --chrg CHARGE
```

The per-conformer xTB property command is equivalent to:

```text
xtb --gbsa toluene --lmo --vfukui --esp -P THREADS --chrg CHARGE input.xyz
```

The backend passes the pinned xTB executable to CREST explicitly. `--quick`
is available only as a smoke-test profile and is recorded as non-production
in the provenance artifact.

## LMO centre convention

For each phosphorus donor, StericX:

1. reads `lmocent.coord` and converts Bohr coordinates to Ångströms;
2. identifies the three nearest non-hydrogen phosphorus substituents;
3. considers the four LMO centres nearest phosphorus;
4. chooses the centre whose minimum distance to those substituents is largest;
5. normalizes the phosphorus-to-LMO vector to a 2.1 Å virtual-metal distance.

This reproduces the free-phosphine LMO selection in Kraken rather than
inferring the coordination direction from molecular geometry.

## Cache and provenance contract

Cache keys include input coordinates, calculation settings, tool versions, and
executable hashes. CREST ensemble generation and xTB LMO properties are
separate stages:

```text
input XYZ
   └── crest/<content key>/       durable conformers + populations
          ├── lmo/<content key>/  independent property result per geometry
          └── jobs/<content key>/ resumable join state + final ensemble
```

The CREST stage is atomically promoted immediately after normal termination
and successful parsing, before any LMO calculation starts. LMO work is
deduplicated by coordinate hash and dispatched through a worker pool bounded
by the requested process count, threads per xTB process, and available CPUs.
The job state is atomically checkpointed after each LMO completion, so an
interrupted retry reuses both the complete CREST ensemble and every completed
LMO cache entry.

Lock files record a UUID, hostname, PID, Linux process-start token, thread, and
creation time. A live local lock cannot be stolen. Locks whose PID is dead or
whose PID start token changed are reclaimed; foreign-host locks require the
configured age threshold. Owners verify both inode and UUID during release so
they cannot remove a replacement lock.

Completed cache entries are immutable. Each manifest records the exact
argument vector, environment controls, input/output hashes, runtime, and logs.
Scratch calculation directories live below `.stericx/cache` and are atomically
promoted only after successful validation. The schema-v1 combined smoke cache
is migrated into the split CREST stage with its original manifest hash and
artifact provenance preserved.
