from constants import G, SIDEREAL_YEAR_TO_SIDEREAL_DAY


class CelestialBody:
    def __init__(self, name, mass=None, radius=None, oblateness=None, J2=None, SOI=None,
                 semimajor=None, T_sidereal=None):

        # given quantities
        self.name = name  # `str` (--) name of the body
        self.mass = mass  # `float` (kg) mass
        self.radius_equatorial = radius  # `float` (km) equatorial radius
        self.oblateness = oblateness  # `float` (--) oblateness
        self.J2 = J2  # `float` (--) second harmonic (J2) parameter
        self.SOI = SOI  # `float` (km) radius of Sphere of Influence (SOI)
        
        # TODO: implement this and alias below
        # self.rotational_rate = rotational_rate  # `float` (rad/s) rotational rate about axis
        
        self.primary = None  # `CelestialBody` (--) the body's primary (if it has one)
        self.semimajor = semimajor  # `float` (km) orbital semi-major axis about primary
        self.T_siderial = T_sidereal  # `float` (day) orbital sidereal period about primary
        
        self.satellites = []  # `list` of `CelestialBody` (--) the satellites orbiting this body

        self.mu = G*self.mass  # `float` (km**3/s**2) gravitational parameter of the body

        # aliases
        self.M = self.mass
        self.R = self.radius_equatorial
        # self.w = self.rotational_rates
        self.a = self.semimajor
        self.T = self.T_siderial


''' Data for major celestial BODIES in the Solar System. Data from 
Orbital Mechanics for Engineering Students, 4 ed, Curtis. Includes
data from Tables 4.3, A.1. '''
BODIES = {
    'Sun'    : CelestialBody('Sun',     1.989e30, 696e3),
    'Mercury': CelestialBody('Mercury', 330.0e21, 2440,  0.0,      60e-6,      112e3,  57.91e6, 87.97),
    'Venus'  : CelestialBody('Venus',   4.869e24, 6052,  0.0,      4.458e-6,   616e3,  108.2e6, 224.7),
    'Earth'  : CelestialBody('Earth',   5.974e24, 6378,  0.003353, 1.08263e-3, 925e3,  149.6e6, 365.256),
    'Moon'   : CelestialBody('Moon',    73.48e21, 1737,  0.0012,   202.7e-6,   66100,  384.4e3, 27.322),
    'Mars'   : CelestialBody('Mars',    641.9e21, 3396,  0.00648,  1.96045e-3, 577e3,  227.9e6, 1.881*SIDEREAL_YEAR_TO_SIDEREAL_DAY),
    'Jupiter': CelestialBody('Jupiter', 1.899e27, 71490, 0.06487,  14.736e-3,  48.2e6, 778.6e6, 11.86*SIDEREAL_YEAR_TO_SIDEREAL_DAY),
    'Saturn' : CelestialBody('Saturn',  568.5e24, 60270, 0.09796,  16.298e-3,  54.8e6, 1.433e9, 29.46*SIDEREAL_YEAR_TO_SIDEREAL_DAY),
    'Uranus' : CelestialBody('Uranus',  86.83e24, 25560, 0.02293,  3.34343e-3, 54.8e6, 2.872e9, 84.01*SIDEREAL_YEAR_TO_SIDEREAL_DAY),
    'Neptune': CelestialBody('Neptune', 102.4e24, 24764, 0.0012,   202.7e-6,   86.6e6, 4.495e9, 164.8*SIDEREAL_YEAR_TO_SIDEREAL_DAY),
    'Pluto'  : CelestialBody('Pluto',   13.13e21, 1187, SOI=3.08e6, semimajor=5.906e9, T_sidereal=247.9*SIDEREAL_YEAR_TO_SIDEREAL_DAY)
}

# set primaries and satellites
PLANET_NAMES = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
for planet in PLANET_NAMES:
    BODIES[planet].primary = BODIES['Sun']
    BODIES['Sun'].satellites.append(BODIES[planet])
BODIES['Earth'].satellites.append(BODIES['Moon'])
BODIES['Moon'].primary = BODIES['Earth']
