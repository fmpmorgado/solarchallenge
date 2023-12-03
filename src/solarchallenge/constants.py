"""
Module to store constants
"""
from datetime import datetime
import numpy as np

#Boltzmann constant in eV/K
BOLTZMANN = 8.617332478E-05

#Stefan-Boltzmann constant [W / (m**2 K**4)]
STEFAN_BOLTZMANN = 5.670374419E-8 

#J2000 time [datetime] in UTC format
J2000 = datetime(2000, 1, 1, 11, 58, 55, 816)

#Earth time of perihelion [datetime in UTC]
T_PERI_EARTH = datetime(2000,1,4,16,17)

#Earth standard gravitational parameter [km**3/s**2]
MU_EARTH = 3.986E5 

#Mean Radius of Earth [km]
R_MEAN_EARTH = 6_371

#Mean albedo factor of Earth
ALBEDO_EARTH = 0.3 

#Mean Earth IR emission [W/m**2]
IR_EARTH = 237.0

#Sun standard gravitational parameter [km**3/s**2]
MU_SUN = 1.327E11

#Mean Sun radius [km]
R_MEAN_SUN = 696_340

#Sun surface area [km ** 2]
A_SURF_SUN = 4 * np.pi * R_MEAN_SUN ** 2

#Sun surface temperature [K]
T_SURF_SUN = 5_778
