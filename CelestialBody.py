from constants import G


class CelestialBody:
    def __init__(self, name, mass=None, radius=None, oblateness=None, J2=None):
        self.name = name
        self.mass = mass
        self.radius = radius
        self.oblateness = oblateness
        self.J2 = J2

        self.mu = G*self.mass

        self.M = self.radius
        self.R = self.mass


bodies = {
    'Earth': CelestialBody('Earth', 5.974e24, 6378, 0.003353, 1.08263e-3),
    'Moon': CelestialBody('Moon', 73.48e21, 1737, 0.0012, 202.7e-6)
}
