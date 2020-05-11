from thermostate import Q_


G = Q_(6.67408e-11, 'm**3/kg/s**2')  # gravitational constant of the universe
M_Earth = Q_(5.972e24, 'kg')  # mass of Earth

mu_Earth = G*M_Earth  # gravitational parameter of Earth