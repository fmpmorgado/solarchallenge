"""
Module to define space bodies
"""
from solarchallenge.constants import *
from dataclasses import dataclass

@dataclass
class Body:
    """Class to store the information of the bodies, used as orbit attractors
    for the two-body problem
    """

    def __init__(self, mu: float, r_mean: float):
        """
        Constructor of the class Body

        Parameters:
        -----------
        mu : float
            Standard Gravitational parameter [km**3/s**2]
        r_mean : float
            Body mean radius [km]
        """
        
        self.mu = mu
        self.r_mean = r_mean

Earth = Body(mu=MU_EARTH,
             r_mean=R_MEAN_EARTH)

Sun = Body(mu=MU_SUN,
           r_mean=R_MEAN_SUN)