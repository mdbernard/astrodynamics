import numpy as np
from constants import G


class CelestialBody:
    def __init__(self, name, mass=None, radius=None, oblateness=None, J2=None, SOI=None,
                 rotational_rate=None, primary=None, satellites=[]):

        # given quantities
        self.name = name  # `str` (--) name of the body
        self.mass = mass  # `float` (kg) mass
        self.radius_equatorial = radius  # `float` (km) equatorial radius
        self.oblateness = oblateness  # `float` (--) oblateness
        self.J2 = J2  # `float` (--) second harmonic (J2) parameter
        self.SOI = SOI  # `float` (km) radius of Sphere of Influence (SOI)
        self.rotational_rate = rotational_rate  # `float` (rad/s) rotational rate about axis
        self.primary = primary  # `str` (--) name of the body's primary (if it has one)
        self.satellites = satellites  # `list` of `CelestialBody` (--) the bodies orbiting this body

        self.mu = G*self.mass  # `float` (km**3/s**2) gravitational parameter of the body

        # aliases
        self.M = self.mass
        self.R = self.radius_equatorial
        self.w = None


''' Data for major celestial bodies in the Solar System. Data from 
Orbital Mechanics for Engineering Students, 4 ed, Curtis. Includes
data from Tables 4.3, A.1.
'''
bodies = {
    'Sun': CelestialBody('Sun', 1.989e30, 696e3),
    'Mercury': CelestialBody('Mercury', 330.0e21, 2440, 0.0, 60e-6, 112e3),
    'Venus': CelestialBody('Venus', 4.869e24, 6052, 0.0, 4.458e-6, 616e3),
    'Earth': CelestialBody('Earth', 5.974e24, 6378, 0.003353, 1.08263e-3, 7.292115e-5, 925e3),
    'Moon': CelestialBody('Moon', 73.48e21, 1737, 0.0012, 202.7e-6, 66100),
    'Mars': CelestialBody('Mars', 641.9e21, 3396, 0.00648, 1.96045e-3, 577e3),
    'Jupiter': CelestialBody('Jupiter', 1.899e27, 71490, 0.06487, 14.736e-3, 48.2e6),
    'Saturn': CelestialBody('Saturn', 568.5e24, 60270, 0.09796, 16.298e-3, 54.8e6),
    'Uranus': CelestialBody('Uranus', 86.83e24, 25560, 0.02293, 3.34343e-3, 54.8e6),
    'Neptune': CelestialBody('Neptune', 102.4e24, 24764, 0.0012, 202.7e-6, 86.6e6),
    'Pluto': CelestialBody('Pluto', 13.13e21, 1187, SOI=3.08e6)
}

# set primaries and satellites
planets = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
for planet in planets:
    bodies[planet].primary = bodies['Sun']
    bodies['Sun'].satellites.append(bodies[planet])
bodies['Earth'].satellites.append(bodies['Moon'])
bodies['Moon'].primary = bodies['Earth']
