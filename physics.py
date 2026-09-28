import math

G = 6.67430e-11
M_EARTH = 5.97219e24
R_EARTH = 6_371_000.0
MU_EARTH = G * M_EARTH


def norm(x, y):
    return math.hypot(x, y)


def gravity_acceleration(x, y, mu=MU_EARTH):
    r = norm(x, y)
    if r == 0:
        raise ValueError("Position cannot be at Earth's center.")
    factor = -mu / (r ** 3)
    return factor * x, factor * y


def exponential_density(altitude_m, rho0=1.225, scale_height_m=8500.0):
    if altitude_m < 0:
        altitude_m = 0
    return rho0 * math.exp(-altitude_m / scale_height_m)


def drag_acceleration(vx, vy, altitude_m, mass_kg, cd=0.47, area_m2=0.01):
    if mass_kg <= 0:
        raise ValueError("mass_kg must be positive")
    speed = norm(vx, vy)
    if speed == 0:
        return 0.0, 0.0
    rho = exponential_density(altitude_m)
    k = 0.5 * rho * cd * area_m2 / mass_kg
    return -k * speed * vx, -k * speed * vy


def derivatives(state, use_drag=False, mass_kg=1.0, cd=0.47, area_m2=0.01):
    x, y, vx, vy = state
    ax_g, ay_g = gravity_acceleration(x, y)
    ax, ay = ax_g, ay_g
    if use_drag:
        altitude = norm(x, y) - R_EARTH
        ax_d, ay_d = drag_acceleration(vx, vy, altitude, mass_kg, cd, area_m2)
        ax += ax_d
        ay += ay_d
    return (vx, vy, ax, ay)


def euler_step(state, dt, deriv_func, **kwargs):
    d = deriv_func(state, **kwargs)
    return tuple(s + dt * ds for s, ds in zip(state, d))


def rk4_step(state, dt, deriv_func, **kwargs):
    k1 = deriv_func(state, **kwargs)
    s2 = tuple(s + 0.5 * dt * k for s, k in zip(state, k1))
    k2 = deriv_func(s2, **kwargs)
    s3 = tuple(s + 0.5 * dt * k for s, k in zip(state, k2))
    k3 = deriv_func(s3, **kwargs)
    s4 = tuple(s + dt * k for s, k in zip(state, k3))
    k4 = deriv_func(s4, **kwargs)
    return tuple(
        s + (dt / 6.0) * (a + 2*b + 2*c + d)
        for s, a, b, c, d in zip(state, k1, k2, k3, k4)
    )


def specific_energy(state, mu=MU_EARTH):
    x, y, vx, vy = state
    r = norm(x, y)
    return 0.5 * (vx*vx + vy*vy) - mu / r
