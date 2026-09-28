# ASTRA v1.0

**Computational Astrodynamics & Trajectory Analysis Laboratory**

ASTRA is a Python computational-physics project for studying trajectories, orbital mechanics, perturbations, numerical integration, and sensitivity to initial conditions.

## What ASTRA can do

- Simulate 2D Earth orbits with Newtonian gravity.
- Compare **Euler** and **RK4** integration and measure energy drift.
- Model a simple exponential atmosphere with aerodynamic drag.
- Add an equatorial-plane **J2 perturbation** term.
- Simulate projectiles with optional **Earth-rotation launch boost**.
- Run a normalized **Earth-Moon circular restricted three-body problem (CR3BP)**.
- Perform **sensitivity analysis** by perturbing initial orbital velocity.
- Export reproducible CSV data and SVG trajectory plots.
- Run automated unit tests and an integrator benchmark experiment.

## Why the project matters

Numerical simulations approximate continuous physics using finite time steps. ASTRA shows how a numerical method can change a trajectory even when the underlying physical equations are identical. This makes the project both an astrodynamics simulator and a numerical-analysis laboratory.

## Quick start

```bash
python astra.py
```

This runs a 400 km orbit with both Euler and RK4.

### Orbit + J2

```bash
python astra.py --mode orbit --altitude-km 400 --dt 10 --j2
```

### Projectile + drag + Earth rotation

```bash
python astra.py --mode projectile --speed-mps 800 --angle-deg 45 --dt 0.1 --drag --earth-rotation --latitude-deg 24.45
```

### Earth-Moon restricted three-body model

```bash
python astra.py --mode threebody --dt 0.001 --duration 10 --method rk4
```

The CR3BP mode uses standard normalized rotating-frame units.

### Sensitivity experiment

```bash
python astra.py --mode sensitivity --altitude-km 400 --dt 10 --perturbation 0.0001
```

## Run tests

```bash
python -m unittest discover -s tests -v
```

## Run benchmark experiment

```bash
python experiments.py
```

This produces `experiments/results/integrator_benchmark.csv` and `benchmark.md`.

## Project structure

```text
ASTRA/
├── astra.py              # CLI + simulations + CSV/SVG output
├── physics.py            # physics models and numerical integrators
├── experiments.py        # numerical-method benchmark
├── tests/
│   └── test_astra.py     # automated tests
├── experiments/
│   └── results/          # generated benchmark reports
├── output/               # generated simulation results
├── LICENSE
└── README.md
```

## Physics implemented

### Two-body gravity

`a = -mu * r_vec / r^3`

### Atmospheric drag

`a_D = -(1/2) * rho * C_D * (A/m) * v * v_vec`

with a simplified exponential atmosphere.

### J2 perturbation

ASTRA includes a 2D equatorial approximation to Earth's oblateness perturbation. It is intentionally a teaching model, not a mission-grade propagator.

### CR3BP

The planar Earth-Moon circular restricted three-body problem is implemented in normalized rotating coordinates, allowing qualitative investigation of multi-body dynamics.

## Numerical methods

- Forward Euler — simple, fast, but accumulates large orbital error.
- Classical RK4 — more computational work per step, but dramatically better short-term orbital accuracy.

## Limitations

ASTRA v1.0 is educational research software, not a spacecraft navigation or mission-design tool. The atmospheric model is simplified, J2 is restricted to a 2D approximation, and CR3BP uses idealized circular normalized dynamics.

## Author

**Hasan Alhashmi (Astro7asan)** — 2026

## License

MIT License.
