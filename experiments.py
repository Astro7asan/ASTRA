import csv
import math
import time
from pathlib import Path

from astra import simulate, orbit_state, energy_drift
from physics import MU_EARTH


def benchmark(output_dir="experiments/results"):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    initial = orbit_state(400)
    period = 2 * math.pi * math.sqrt(initial[0]**3 / MU_EARTH)
    rows = []
    for method in ("euler", "rk4"):
        for dt in (1, 5, 10, 30, 60):
            t0 = time.perf_counter()
            sim = simulate(initial, dt, period, method=method)
            elapsed = time.perf_counter() - t0
            rows.append({
                "method": method,
                "dt_s": dt,
                "samples": len(sim),
                "energy_drift_percent": energy_drift(sim),
                "runtime_s": elapsed,
                "final_altitude_km": sim[-1]["altitude_m"] / 1000,
            })

    with (out/"integrator_benchmark.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    md = [
        "# ASTRA Integrator Benchmark",
        "",
        "| Method | dt (s) | Samples | Energy drift (%) | Runtime (s) | Final altitude (km) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        md.append(
            f"| {r['method']} | {r['dt_s']} | {r['samples']} | "
            f"{r['energy_drift_percent']:.6g} | {r['runtime_s']:.6f} | "
            f"{r['final_altitude_km']:.3f} |"
        )
    (out/"benchmark.md").write_text("\n".join(md), encoding="utf-8")
    return rows


if __name__ == "__main__":
    benchmark()
