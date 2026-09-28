import argparse
import csv
import math
from pathlib import Path

from physics import R_EARTH, MU_EARTH, derivatives, euler_step, rk4_step, norm, specific_energy


def simulate(initial_state, dt, duration, method="rk4", use_drag=False,
             mass_kg=1.0, cd=0.47, area_m2=0.01, stop_at_ground=False):
    stepper = rk4_step if method == "rk4" else euler_step
    state = tuple(float(v) for v in initial_state)
    t = 0.0
    rows = []
    while t <= duration + 1e-12:
        x, y, vx, vy = state
        r = norm(x, y)
        rows.append({
            "t_s": t,
            "x_m": x,
            "y_m": y,
            "vx_mps": vx,
            "vy_mps": vy,
            "altitude_m": r - R_EARTH,
            "speed_mps": norm(vx, vy),
            "specific_energy_J_per_kg": specific_energy(state),
        })
        if stop_at_ground and t > 0 and r <= R_EARTH:
            break
        state = stepper(
            state, dt, derivatives,
            use_drag=use_drag, mass_kg=mass_kg, cd=cd, area_m2=area_m2
        )
        t += dt
    return rows


def save_csv(rows, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def save_svg(series, path, width=900, height=700, margin=50):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    all_points = [(row["x_m"], row["y_m"]) for _, rows in series for row in rows]
    xs = [p[0] for p in all_points]
    ys = [p[1] for p in all_points]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    span_x = max(max_x - min_x, 1.0)
    span_y = max(max_y - min_y, 1.0)
    scale = min((width - 2*margin)/span_x, (height - 2*margin)/span_y)

    def tx(x): return margin + (x - min_x) * scale
    def ty(y): return height - margin - (y - min_y) * scale

    palette = ["#3b82f6", "#ef4444", "#10b981", "#f59e0b"]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="20" y="30" font-family="sans-serif" font-size="20">ASTRA Trajectory Comparison</text>'
    ]
    for i, (label, rows) in enumerate(series):
        pts = " ".join(f"{tx(r['x_m']):.2f},{ty(r['y_m']):.2f}" for r in rows)
        color = palette[i % len(palette)]
        lines.append(f'<polyline fill="none" stroke="{color}" stroke-width="2" points="{pts}"/>')
        lines.append(f'<text x="20" y="{60 + 22*i}" font-family="sans-serif" font-size="14" fill="{color}">{label}</text>')
    lines.append('</svg>')
    path.write_text("\n".join(lines), encoding="utf-8")


def circular_orbit_state(altitude_km):
    r = R_EARTH + altitude_km * 1000.0
    v = math.sqrt(MU_EARTH / r)
    return (r, 0.0, 0.0, v)


def projectile_state(altitude_km, speed_mps, angle_deg):
    r = R_EARTH + altitude_km * 1000.0
    a = math.radians(angle_deg)
    # Local vertical is +x, local horizontal is +y at the launch point.
    vx = speed_mps * math.sin(a)
    vy = speed_mps * math.cos(a)
    return (r, 0.0, vx, vy)


def parse_args():
    p = argparse.ArgumentParser(description="ASTRA: 2D astrodynamics and trajectory simulator")
    p.add_argument("--mode", choices=["orbit", "projectile"], default="orbit")
    p.add_argument("--altitude-km", type=float, default=400.0)
    p.add_argument("--dt", type=float, default=10.0)
    p.add_argument("--duration", type=float, default=None)
    p.add_argument("--method", choices=["euler", "rk4", "both"], default="both")
    p.add_argument("--speed-mps", type=float, default=1000.0, help="Projectile launch speed")
    p.add_argument("--angle-deg", type=float, default=45.0, help="Projectile angle above local horizontal")
    p.add_argument("--drag", action="store_true")
    p.add_argument("--mass-kg", type=float, default=10.0)
    p.add_argument("--cd", type=float, default=0.47)
    p.add_argument("--area-m2", type=float, default=0.01)
    p.add_argument("--output-dir", default="sample_output")
    return p.parse_args()


def main():
    args = parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    if args.mode == "orbit":
        initial = circular_orbit_state(args.altitude_km)
        r = initial[0]
        period = 2 * math.pi * math.sqrt(r**3 / MU_EARTH)
        duration = args.duration if args.duration is not None else period
        stop_at_ground = False
    else:
        initial = projectile_state(args.altitude_km, args.speed_mps, args.angle_deg)
        duration = args.duration if args.duration is not None else 300.0
        stop_at_ground = True

    methods = ["euler", "rk4"] if args.method == "both" else [args.method]
    series = []
    for method in methods:
        rows = simulate(
            initial, args.dt, duration, method=method,
            use_drag=args.drag, mass_kg=args.mass_kg, cd=args.cd, area_m2=args.area_m2,
            stop_at_ground=stop_at_ground,
        )
        save_csv(rows, out / f"{args.mode}_{method}.csv")
        series.append((method.upper(), rows))
        e0 = rows[0]["specific_energy_J_per_kg"]
        e1 = rows[-1]["specific_energy_J_per_kg"]
        drift = abs((e1 - e0) / e0) * 100 if e0 else float("nan")
        print(f"{method.upper()}: {len(rows)} samples, final altitude={rows[-1]['altitude_m']/1000:.3f} km, energy drift={drift:.6f}%")

    save_svg(series, out / "trajectory.svg")
    print(f"Saved results to: {out.resolve()}")


if __name__ == "__main__":
    main()
