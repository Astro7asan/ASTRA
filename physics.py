import math

G = 6.67430e-11
M_EARTH = 5.97219e24
R_EARTH = 6_371_000.0
MU_EARTH = G * M_EARTH
J2 = 1.08262668e-3
M_MOON = 7.342e22
MU_MOON = G * M_MOON
EARTH_MOON_DISTANCE = 384_400_000.0
EARTH_ROTATION_RATE = 7.2921159e-5


def norm(*values):
    return math.sqrt(sum(v * v for v in values))


def gravity_acceleration(x, y, mu=MU_EARTH):
    r = norm(x, y)
    if r == 0:
        raise ValueError("Position cannot be at the gravitating body's center.")
    f = -mu / (r ** 3)
    return f * x, f * y


def j2_acceleration_2d(x, y):
    """Equatorial-plane J2 correction for a 2D Earth-centered model."""
    r = norm(x, y)
    if r == 0:
        raise ValueError("Position cannot be at Earth's center.")
    factor = -1.5 * J2 * MU_EARTH * (R_EARTH ** 2) / (r ** 5)
    return factor * x, factor * y


def exponential_density(altitude_m, rho0=1.225, scale_height_m=8500.0):
    altitude_m = max(0.0, altitude_m)
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


def derivatives(state, use_drag=False, use_j2=False, mass_kg=1.0, cd=0.47, area_m2=0.01):
    x, y, vx, vy = state
    ax, ay = gravity_acceleration(x, y)
    if use_j2:
        jx, jy = j2_acceleration_2d(x, y)
        ax += jx
        ay += jy
    if use_drag:
        altitude = norm(x, y) - R_EARTH
        dx, dy = drag_acceleration(vx, vy, altitude, mass_kg, cd, area_m2)
        ax += dx
        ay += dy
    return vx, vy, ax, ay


def cr3bp_derivatives(state, mu_ratio):
    """Planar circular restricted three-body problem in normalized rotating units."""
    x, y, vx, vy = state
    mu = mu_ratio
    r1 = math.sqrt((x + mu) ** 2 + y ** 2)
    r2 = math.sqrt((x - 1 + mu) ** 2 + y ** 2)
    if r1 == 0 or r2 == 0:
        raise ValueError("State intersects a primary body.")
    ax = x + 2 * vy - (1 - mu) * (x + mu) / (r1 ** 3) - mu * (x - 1 + mu) / (r2 ** 3)
    ay = y - 2 * vx - (1 - mu) * y / (r1 ** 3) - mu * y / (r2 ** 3)
    return vx, vy, ax, ay


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
    return tuple(s + (dt / 6.0) * (a + 2*b + 2*c + d) for s, a, b, c, d in zip(state, k1, k2, k3, k4))


def specific_energy(state, mu=MU_EARTH):
    x, y, vx, vy = state
    r = norm(x, y)
    return 0.5 * (vx * vx + vy * vy) - mu / r


def angular_momentum_z(state):
    x, y, vx, vy = state
    return x * vy - y * vx


def orbital_elements_2d(state, mu=MU_EARTH):
    x, y, vx, vy = state
    r = norm(x, y)
    v2 = vx * vx + vy * vy
    h = angular_momentum_z(state)
    energy = 0.5 * v2 - mu / r
    a = math.inf if abs(energy) < 1e-20 else -mu / (2 * energy)
    e = math.sqrt(max(0.0, 1 + 2 * energy * h * h / (mu * mu)))
    return {"semi_major_axis_m": a, "eccentricity": e, "specific_energy": energy, "h_z": h}


def circular_orbit_speed(altitude_m):
    return math.sqrt(MU_EARTH / (R_EARTH + altitude_m))


def earth_rotation_surface_speed(latitude_deg):
    return EARTH_ROTATION_RATE * R_EARTH * math.cos(math.radians(latitude_deg))
