import argparse
import csv
import math
from pathlib import Path

from physics import (
    R_EARTH, MU_EARTH, M_EARTH, M_MOON, derivatives, cr3bp_derivatives,
    euler_step, rk4_step, norm, specific_energy, orbital_elements_2d,
    circular_orbit_speed, earth_rotation_surface_speed,
)


def simulate(initial_state, dt, duration, method="rk4", deriv_func=derivatives,
             stop_at_ground=False, deriv_kwargs=None):
    deriv_kwargs = deriv_kwargs or {}
    stepper = rk4_step if method == "rk4" else euler_step
    state = tuple(float(v) for v in initial_state)
    rows, t = [], 0.0
    while t <= duration + 1e-12:
        x, y, vx, vy = state
        r = norm(x, y)
        row = {"t_s": t, "x": x, "y": y, "vx": vx, "vy": vy, "r": r, "speed": norm(vx, vy)}
        if deriv_func is derivatives:
            row["altitude_m"] = r - R_EARTH
            row["specific_energy_J_per_kg"] = specific_energy(state)
        rows.append(row)
        if stop_at_ground and t > 0 and r <= R_EARTH:
            break
        state = stepper(state, dt, deriv_func, **deriv_kwargs)
        t += dt
    return rows


def save_csv(rows, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def save_svg(series, path, title="ASTRA trajectory", width=900, height=700):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pts_all = [(r["x"], r["y"]) for _, rows in series for r in rows]
    xs, ys = [p[0] for p in pts_all], [p[1] for p in pts_all]
    min_x, max_x, min_y, max_y = min(xs), max(xs), min(ys), max(ys)
    sx, sy = max(max_x-min_x, 1e-12), max(max_y-min_y, 1e-12)
    margin = 60
    scale = min((width-2*margin)/sx, (height-2*margin)/sy)
    tx = lambda x: margin + (x-min_x)*scale
    ty = lambda y: height - margin - (y-min_y)*scale
    colors = ["#2563eb", "#dc2626", "#059669", "#d97706", "#7c3aed"]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="20" y="30" font-family="sans-serif" font-size="20">{title}</text>'
    ]
    for i, (label, rows) in enumerate(series):
        pts = " ".join(f"{tx(r['x']):.2f},{ty(r['y']):.2f}" for r in rows)
        c = colors[i % len(colors)]
        lines.append(f'<polyline fill="none" stroke="{c}" stroke-width="2" points="{pts}"/>')
        lines.append(f'<text x="20" y="{60+20*i}" font-family="sans-serif" font-size="13" fill="{c}">{label}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def orbit_state(altitude_km, speed_scale=1.0):
    r = R_EARTH + altitude_km * 1000
    return (r, 0.0, 0.0, circular_orbit_speed(altitude_km * 1000) * speed_scale)


def projectile_state(altitude_km, speed, angle_deg, latitude_deg=0.0, include_rotation=False):
    r = R_EARTH + altitude_km * 1000
    a = math.radians(angle_deg)
    vx = speed * math.sin(a)
    vy = speed * math.cos(a)
    if include_rotation:
        vy += earth_rotation_surface_speed(latitude_deg)
    return (r, 0.0, vx, vy)


def energy_drift(rows):
    e0 = rows[0]["specific_energy_J_per_kg"]
    e1 = rows[-1]["specific_energy_J_per_kg"]
    return abs((e1-e0)/e0) * 100 if e0 else float("nan")


def run_orbit(args):
    initial = orbit_state(args.altitude_km, args.speed_scale)
    period = 2 * math.pi * math.sqrt(initial[0]**3 / MU_EARTH)
    duration = args.duration if args.duration is not None else period
    methods = ["euler", "rk4"] if args.method == "both" else [args.method]
    series = []
    for method in methods:
        rows = simulate(
            initial, args.dt, duration, method=method,
            deriv_kwargs={
                "use_drag": args.drag, "use_j2": args.j2,
                "mass_kg": args.mass_kg, "cd": args.cd, "area_m2": args.area_m2
            }
        )
        save_csv(rows, Path(args.output_dir)/f"orbit_{method}.csv")
        series.append((method.upper(), rows))
        final_state = (rows[-1]["x"], rows[-1]["y"], rows[-1]["vx"], rows[-1]["vy"])
        elems = orbital_elements_2d(final_state)
        print(f"{method.upper()}: final altitude={rows[-1]['altitude_m']/1000:.3f} km | energy drift={energy_drift(rows):.6g}% | e={elems['eccentricity']:.6f}")
    save_svg(series, Path(args.output_dir)/"trajectory.svg", "ASTRA orbit comparison")


def run_projectile(args):
    initial = projectile_state(args.altitude_km, args.speed_mps, args.angle_deg, args.latitude_deg, args.earth_rotation)
    duration = args.duration or 300.0
    rows = simulate(
        initial, args.dt, duration,
        method=args.method if args.method != "both" else "rk4",
        stop_at_ground=True,
        deriv_kwargs={
            "use_drag": args.drag, "use_j2": False,
            "mass_kg": args.mass_kg, "cd": args.cd, "area_m2": args.area_m2
        }
    )
    save_csv(rows, Path(args.output_dir)/"projectile.csv")
    save_svg([("PROJECTILE", rows)], Path(args.output_dir)/"trajectory.svg", "ASTRA projectile")
    print(f"Impact/finish time={rows[-1]['t_s']:.2f} s | max altitude={max(r['altitude_m'] for r in rows)/1000:.3f} km")


def run_threebody(args):
    mu = M_MOON / (M_EARTH + M_MOON)
    initial = (args.x0, args.y0, args.vx0, args.vy0)
    rows = simulate(
        initial, args.dt, args.duration or 10.0,
        method=args.method if args.method != "both" else "rk4",
        deriv_func=cr3bp_derivatives,
        deriv_kwargs={"mu_ratio": mu}
    )
    save_csv(rows, Path(args.output_dir)/"threebody.csv")
    save_svg([("CR3BP", rows)], Path(args.output_dir)/"trajectory.svg", "ASTRA Earth-Moon CR3BP")
    print(f"Earth-Moon CR3BP completed: {len(rows)} samples")


def run_sensitivity(args):
    base = orbit_state(args.altitude_km, 1.0)
    pert = orbit_state(args.altitude_km, 1.0 + args.perturbation)
    duration = args.duration or 2*math.pi*math.sqrt(base[0]**3/MU_EARTH)
    a = simulate(base, args.dt, duration, method="rk4")
    b = simulate(pert, args.dt, duration, method="rk4")
    save_csv(a, Path(args.output_dir)/"sensitivity_base.csv")
    save_csv(b, Path(args.output_dir)/"sensitivity_perturbed.csv")
    sep = norm(a[-1]["x"]-b[-1]["x"], a[-1]["y"]-b[-1]["y"])
    save_svg([("BASE", a), ("PERTURBED", b)], Path(args.output_dir)/"trajectory.svg", "ASTRA sensitivity analysis")
    print(f"Initial velocity perturbation={args.perturbation*100:.6f}% | final position separation={sep/1000:.3f} km")


def parser():
    p = argparse.ArgumentParser(description="ASTRA v1.0 - Computational Astrodynamics & Trajectory Analysis Laboratory")
    p.add_argument("--mode", choices=["orbit", "projectile", "threebody", "sensitivity"], default="orbit")
    p.add_argument("--method", choices=["euler", "rk4", "both"], default="both")
    p.add_argument("--dt", type=float, default=10.0)
    p.add_argument("--duration", type=float)
    p.add_argument("--output-dir", default="output")
    p.add_argument("--altitude-km", type=float, default=400.0)
    p.add_argument("--speed-scale", type=float, default=1.0)
    p.add_argument("--drag", action="store_true")
    p.add_argument("--j2", action="store_true")
    p.add_argument("--mass-kg", type=float, default=10.0)
    p.add_argument("--cd", type=float, default=0.47)
    p.add_argument("--area-m2", type=float, default=0.01)
    p.add_argument("--speed-mps", type=float, default=800.0)
    p.add_argument("--angle-deg", type=float, default=45.0)
    p.add_argument("--latitude-deg", type=float, default=24.45)
    p.add_argument("--earth-rotation", action="store_true")
    p.add_argument("--x0", type=float, default=0.8)
    p.add_argument("--y0", type=float, default=0.0)
    p.add_argument("--vx0", type=float, default=0.0)
    p.add_argument("--vy0", type=float, default=0.25)
    p.add_argument("--perturbation", type=float, default=1e-4)
    return p


def main():
    args = parser().parse_args()
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    if args.dt <= 0:
        raise SystemExit("dt must be positive")
    {"orbit": run_orbit, "projectile": run_projectile, "threebody": run_threebody, "sensitivity": run_sensitivity}[args.mode](args)
    print(f"Saved outputs to {Path(args.output_dir).resolve()}")


if __name__ == "__main__":
    main()
