# ASTRA Integrator Benchmark

| Method | dt (s) | Samples | Energy drift (%) | Final altitude (km) |
|---|---:|---:|---:|---:|
| Euler | 1 | 5545 | 1.38446 | 495.844 |
| Euler | 5 | 1109 | 6.24051 | 870.343 |
| Euler | 10 | 555 | 11.1594 | 1326.535 |
| Euler | 30 | 185 | 24.1937 | 3099.243 |
| Euler | 60 | 93 | 36.3998 | 5604.251 |
| RK4 | 1 | 5545 | 2.78437e-13 | 400.000 |
| RK4 | 5 | 1109 | 1.01718e-10 | 400.000 |
| RK4 | 10 | 555 | 3.25801e-09 | 400.000 |
| RK4 | 30 | 185 | 7.89019e-07 | 400.000 |
| RK4 | 60 | 93 | 2.52665e-05 | 399.998 |

These results were generated from a 400 km circular-orbit test. Runtime is intentionally omitted from this committed summary because it depends on the computer running the benchmark.
