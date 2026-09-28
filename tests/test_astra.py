import math
import unittest

from astra import simulate, orbit_state
from physics import R_EARTH, MU_EARTH, gravity_acceleration, orbital_elements_2d


class TestASTRA(unittest.TestCase):
    def test_surface_gravity_reasonable(self):
        ax, ay = gravity_acceleration(R_EARTH, 0)
        self.assertAlmostEqual(abs(ax), 9.82, delta=0.05)
        self.assertAlmostEqual(ay, 0.0, delta=1e-12)

    def test_rk4_circular_orbit_stability(self):
        state = orbit_state(400)
        period = 2 * math.pi * math.sqrt(state[0]**3 / MU_EARTH)
        rows = simulate(state, 10, period, method="rk4")
        self.assertLess(abs(rows[-1]["altitude_m"] - 400000), 100.0)

    def test_orbital_elements_near_circular(self):
        e = orbital_elements_2d(orbit_state(400))["eccentricity"]
        self.assertLess(e, 1e-6)

    def test_euler_worse_than_rk4(self):
        state = orbit_state(400)
        period = 2 * math.pi * math.sqrt(state[0]**3 / MU_EARTH)
        euler = simulate(state, 60, period, method="euler")
        rk4 = simulate(state, 60, period, method="rk4")
        de = abs(euler[-1]["specific_energy_J_per_kg"] - euler[0]["specific_energy_J_per_kg"])
        dr = abs(rk4[-1]["specific_energy_J_per_kg"] - rk4[0]["specific_energy_J_per_kg"])
        self.assertGreater(de, dr * 1000)


if __name__ == "__main__":
    unittest.main()
