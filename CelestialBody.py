import numpy as np
from constants import G


class CelestialBody:
    def __init__(self, name, mass=None, radius=None, oblateness=None, J2=None):
        self.name = name
        self.mass = mass
        self.radius_equatorial = radius
        self.oblateness = oblateness
        self.J2 = J2

        self.mu = G*self.mass

        self.M = self.mass
        self.R = self.radius_equatorial

        self.rotational_rate = None
        self.w = None

    def set_rotation_rate(self, rate):
        ''' Set rotational rate about axis. '''
        self.rotational_rate = rate  # (rad/s)
        self.w = self.rotational_rate


bodies = {
    'Earth': CelestialBody('Earth', 5.974e24, 6378, 0.003353, 1.08263e-3),
    'Moon': CelestialBody('Moon', 73.48e21, 1737, 0.0012, 202.7e-6)
}

bodies['Earth'].set_rotation_rate(7.292115e-5)  # (rad/s) Earth rotational velocity
