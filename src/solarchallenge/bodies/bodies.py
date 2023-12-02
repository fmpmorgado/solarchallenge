#https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
from datetime import datetime
class Body:
    def __init__(self, mu): 
        self.mu = mu
#TODO
EARTH = {"mu": 3.986E5,
         "time_perihelion": datetime(2000,1,4,16,17)}

SUN = {"mu": 1.327E11}

Earth = Body(EARTH["mu"])
Sun = Body(SUN["mu"])




"""
Earth = Body(
    k=constants.GM_earth,
    R=constants.R_earth,
    R_mean=constants.R_mean_earth,
    R_polar=constants.R_polar_earth,
    rotational_period=constants.rotational_period_earth,
    mass=constants.M_earth,
    J2=constants.J2_earth,
    J3=constants.J3_earth,
    mean_a=constants.mean_a_earth,
)
"""