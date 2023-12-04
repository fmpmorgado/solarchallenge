"""
Module to store constants
"""
from datetime import datetime, timedelta
import numpy as np


#Astronomical Unit [km]
AU = 149_597_870.7

#Boltzmann constant in eV/K
BOLTZMANN = 8.617332478E-05

#Stefan-Boltzmann constant [W / (m**2 K**4)]
STEFAN_BOLTZMANN = 5.670374419E-8 

#J2000 time [datetime] in UTC format
J2000 = datetime(2000, 1, 1, 11, 58, 55, 816)
S = timedelta(seconds=1)
M = timedelta(minutes=1)
H = timedelta(hours=1)
D = timedelta(days=1)

#Earth time of perihelion [datetime in UTC]
T_PERI_EARTH = datetime(2000, 1, 4, 16, 17)

#Earth standard gravitational parameter [km**3/s**2]
MU_EARTH = 3.986E5 

#Mean Radius of Earth [km]
R_MEAN_EARTH = 6_371

#Mean albedo factor of Earth
ALBEDO_EARTH = 0.3 

#Mean Earth IR emission [W/m**2]
IR_EARTH = 237.0

#Classical orbital Elements of Earth at J2000 (km, deg)
A_EARTH = 1.00000011 * AU
ECC_EARTH = 0.01671022
INC_EARTH = 0.00005
RAAN_EARTH = -11.26064
ARGP_EARTH = 102.94719
NU_EARTH = 356.907

#Sun standard gravitational parameter [km**3/s**2]
MU_SUN = 1.327E11

#Mean Sun radius [km]
R_MEAN_SUN = 696_340

#Sun surface area [km ** 2]
A_SURF_SUN = 4 * np.pi * R_MEAN_SUN ** 2

#Sun surface temperature [K]
T_SURF_SUN = 5_778
