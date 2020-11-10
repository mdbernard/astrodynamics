import numpy as np
from constants import G


class CelestialBody:
    def __init__(self, name, mass=None, radius=None, oblateness=None, J2=None,
                 rotational_rate=None, primary=None, satellites=[]):

        # given quantities
        self.name = name  # `str` (--) name of the body
        self.mass = mass  # `float` (kg) mass
        self.radius_equatorial = radius  # `float` (km) equatorial radius
        self.oblateness = oblateness  # `float` (--) oblateness
        self.J2 = J2  # `float` (--) second harmonic (J2) parameter
        self.rotational_rate = rotational_rate  # `float` (rad/s) rotational rate about axis
        self.primary = primary  # `str` (--) name of the body's primary (if it has one)
        self.satellites = satellites  # `list` of `CelestialBody` (--) the bodies orbiting this body

        self.mu = G*self.mass  # `float` (km**3/s**2) gravitational parameter of the body

        # aliases
        self.M = self.mass
        self.R = self.radius_equatorial
        self.w = None


bodies = {
    'Earth': CelestialBody('Earth', 5.974e24, 6378, 0.003353, 1.08263e-3, 7.292115e-5),
    'Moon': CelestialBody('Moon', 73.48e21, 1737, 0.0012, 202.7e-6)
}
