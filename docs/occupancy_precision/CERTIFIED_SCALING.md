# Certified integration batch scaling

Twelve independent original molecules at certified 0.05 Å³; external Rayon batch pool, no internal parallelism or shared geometry cache. Pool construction/input parsing/one warmup are outside internal batch timing; process metrics include them. Eight alternating worker-count orders plus retained warmup round. Ordered results must equal each frozen certificate exactly.

All 45 processes preserve every individual certificate and output order. Values below are medians of eight measured processes; full ranges, MAD, CPU and RSS distributions are in [certified_scaling.json](certified_scaling.json).

| Workers | Batch ms | MAD ms | Full range ms | Molecules/s | Speedup | Efficiency | Median process peak RSS MiB |
|---|---|---|---|---|---|---|---|
| 1 | 333.128 | 1.396 | 330.786–335.423 | 36.02 | 1.0000× | 100.00% | 16.127 |
| 2 | 171.186 | 0.496 | 170.521–174.081 | 70.10 | 1.9460× | 97.30% | 18.566 |
| 4 | 110.675 | 10.620 | 97.510–125.164 | 108.43 | 3.0100× | 75.25% | 34.650 |
| 6 | 83.392 | 1.682 | 78.110–89.140 | 143.90 | 3.9947× | 66.58% | 42.682 |
| 12 | 89.091 | 0.483 | 84.934–91.226 | 134.69 | 3.7392× | 31.16% | 57.826 |

This is batch-level concurrency for independent molecules. Each certified integration remains serial, with its own candidates, heap and scratch buffers. The pool is created explicitly by the benchmark caller, so no nested pool or new default concurrency policy was added to production. Peak RSS includes the pool and warmup; more simultaneous adaptive heaps can increase memory. Twelve heterogeneous tasks impose a load-balance limit, and the result does not establish scaling of one large integration.
