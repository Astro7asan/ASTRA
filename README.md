# ASTRA

**ASTRA — A Computational Astrodynamics & Trajectory Analysis Laboratory**

ASTRA is a small Python physics project for studying how numerical methods affect simulated trajectories. The starter version models 2D motion around Earth, compares the simple Euler integrator with fourth-order Runge–Kutta (RK4), and can also simulate a projectile with optional atmospheric drag.

## Why this project exists

A computer never follows a continuous trajectory directly. It advances the system in finite time steps, so the numerical method and step size can change the result. ASTRA makes that error visible by running the same physical scenario with different integrators and comparing the resulting paths and energy drift.

## Features

- 2D Newtonian Earth gravity
- Circular-orbit initial conditions
- Projectile mode
- Optional exponential-atmosphere drag
- Euler and RK4 numerical integration
- CSV trajectory export
- SVG trajectory visualization
- Mechanical-energy drift reporting for numerical-method comparison

## Files

- `astra.py` — command-line simulator and output generation
- `physics.py` — constants, force models, Euler, RK4, and energy calculations
- `sample_output/` — generated example results

## Requirements

Python 3.10+ is recommended. This starter version uses only the Python standard library.

## Run it

From the project folder:

```bash
python astra.py
```

Example: 400 km circular orbit with a 60 s time step:

```bash
python astra.py --mode orbit --altitude-km 400 --dt 60
```

Example: projectile with atmospheric drag:

```bash
python astra.py --mode projectile --altitude-km 0.001 --speed-mps 800 --angle-deg 45 --dt 0.1 --duration 120 --drag
```

Choose one integrator:

```bash
python astra.py --method rk4
```

## What to observe

For an orbit, Euler usually accumulates much more numerical error than RK4 at the same time step. Compare the final altitude and reported energy drift in the terminal, then open `sample_output/trajectory.svg` to see the paths.

## Next development steps

The starter project can grow into a larger astrodynamics laboratory with adaptive integration, perturbation models, three-body dynamics, sensitivity analysis, test coverage, notebooks, and an interactive interface.
