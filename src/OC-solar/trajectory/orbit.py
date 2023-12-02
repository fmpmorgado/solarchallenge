import numpy as np
import datetime
from src.bodies.bodies import Body

class Orbit:
    """Class to store the position and velocity of a body at a
    given time (epoch).

    The implicit reference system is an inertial one, based on
    the frame J2000.
    """

    def __init__(self, r: list[float], v: list[float], epoch: datetime.datetime, attractor: Body):
        """Constructor of the class Obit.

        Parameters
        ----------
        r : list[float]
            Position vector [km]
        v : list[float]
            Velocity vector [km]
        epoch : datetime.datetime
            Epoch of the orbit
        attractor : Body
            Object containing information of the attractor (ex: Sun, Earth)
        """

        self.r = r
        self.v = v
        self.epoch = epoch
        self.attractor = attractor

    @classmethod
    def from_vector(cls, r: list[float], v: list[float], epoch: datetime.datetime, attractor: Body):
        """Return Orbit from position and velocity vectors.

        Parameters
        ----------
        r : list[float]
            Position vector [km]
        v : list[float]
            Velocity vector [km]
        epoch : datetime.datetime
            Epoch of the orbit
        attractor : Body
            Object containing information of the attractor (ex: Sun, Earth)
        """

        return cls(r, v, epoch, attractor)

    @classmethod
    def from_coe(cls, a: float, ecc: float, inc: float, raan: float, argp: float, nu: float, epoch: datetime.datetime, attractor: Body):
        """Return Orbit from classical orbital elements.

        Parameters
        ----------
        a : float
            Semi-major axis [km]
        ecc : float 
            Eccentricity
        inc : float
            Inclination [deg]
        raan : float 
            Right ascension of the ascending node [deg]
        argp : float
            Argument of the pericenter [deg]
        nu : float
            True anomaly [deg]
        epoch : Datetime
            Epoch of the orbit
        attractor : Body
            Object containing information of the attractor (ex: Sun, Earth)
        """
        
        r, v = coe2rv(a=a,
                      ecc=ecc,
                      inc=inc,
                      raan=raan,
                      argp=argp,
                      nu=nu,
                      mu=attractor.mu)

        return cls.from_vector(r, v, epoch, attractor)

#https://orbital-mechanics.space/classical-orbital-elements/orbital-elements-and-the-state-vector.html
def coe2rv(a: float, ecc: float, inc: float, raan: float, argp: float, nu: float, mu:float) -> tuple[list[float], list[float]]:
    """Convert classical orbital elements to cartesian

    Parameters
    ----------
    a : float
        Semi-major axis [km]
    ecc : float 
        Eccentricity
    inc : float
        Inclination [deg]
    raan : float 
        Right ascension of the ascending node [deg]
    argp : float
        Argument of the pericenter [deg]
    nu : float
        True anomaly [deg]
    mu : float
        Standard Gravitational Parameter [km**3/s**2]
    
    Returns
    r : list[float]
        Position vector [km]
    v : list[float]
        Velocity vector [km]
    """

    #Convert deg to rad units
    inc *= (np.pi / 180.0)
    raan *= (np.pi / 180.0)
    argp *= (np.pi / 180.0)
    nu *= (np.pi / 180.0)

    #distance from focal point
    p = a * (1 - ecc ** 2)
    d = p / (1 + ecc * np.cos(nu))
    
    #Position in Cartesian frame
    r = [d * (np.cos(raan) * np.cos(argp + nu) - np.sin(raan) * np.sin(argp + nu) * np.cos(inc)),
         d * (np.sin(raan) * np.cos(argp + nu) + np.cos(raan) * np.sin(argp + nu) * np.cos(inc)),
         d * np.sin(argp + nu) * np.sin(inc)]

    #Velocity in Cartesian frame
    v = [- np.sqrt(mu / p) * (np.cos(raan) * (np.sin(argp + nu) + ecc * np.sin(argp)) + 
                              np.sin(raan) * (np.cos(argp + nu) + ecc * np.cos(argp)) * np.cos(inc)),
         - np.sqrt(mu / p) * (np.sin(raan) * (np.sin(argp + nu) + ecc * np.sin(argp)) -
                              np.cos(raan) * (np.cos(argp + nu) + ecc * np.cos(argp)) * np.cos(inc)),
           np.sqrt(mu / p) * (np.cos(argp + nu) + ecc * np.cos(argp)) * np.sin(inc)]

    return r, v

#https://control.asu.edu/Classes/MAE462/462Lecture07.pdf
def rv2coe(r: list[float], v: list[float], mu: float) -> tuple[float, float, float, float, float, float]:
    """Convert cartesian to classical orbital elments

    Parameters
    ----------
    r : list[float]
        Position vector [km]
    v : list[float]
        Velocity vector [km]
    mu : float
        Standard Gravitational Parameter [km**3/s**2]

    Returns
    -------
    a : float
        Semi-major axis [km]
    ecc : float 
        Eccentricity
    inc : float
        Inclination [deg]
    raan : float 
        Right ascension of the ascending node [deg]
    argp : float
        Argument of the pericenter [deg]
    nu : float
        True anomaly [deg]
    """

    # Convert to numpy array
    r = np.array(r)
    v = np.array(v)
    r_norm = np.linalg.norm(r) 
    v_norm = np.linalg.norm(v)

    # Computation of auxiliary variables
    h = np.cross(r,v)
    n = np.cross([0,0,1],h)

    # Eccentricity
    ecc = 1.0 / mu * np.cross(v, h) - r/r_norm
    ecc_norm = np.linalg.norm(ecc) 
    
    # Eccentric anomaly
    e_anom = v_norm ** 2 / 2.0 - mu / r_norm

    # Semi-major axis
    a = -mu / (2 * e_anom)
   
    # Inclination
    inc = np.arccos(np.dot(h / np.linalg.norm(h), [0, 0, 1]))

    # Raan with quadrant correction
    raan = np.arccos(np.dot([1, 0, 0], n / np.linalg.norm(n)))
    if np.dot([0,1,0], n) < 0:
        raan = 2*np.pi - raan

    # Argument of the pericenter with quadrant correction
    argp = np.arccos(np.dot(n, ecc) / (np.linalg.norm(n) * ecc_norm))
    if np.dot([0,0,1], n) < 0:
        argp = 2*np.pi - argp

    # True anomaly with quadrant correction
    nu = np.arccos(np.dot(r, ecc)/(r_norm * ecc_norm))
    if np.dot(r, v) < 0:
        nu = 2*np.pi - nu

    # Convert from deg to rad
    inc /= np.pi/180.0
    raan /= np.pi/180.0
    argp /= np.pi/180.0
    nu /= np.pi/180.0

    return a, ecc_norm, inc, raan, argp, nu