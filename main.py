from src.solarchallenge.trajectory.orbit import Orbit, coe2rv
from src.solarchallenge.bodies.bodies import Earth, Sun
import numpy as np
from datetime import datetime
from solarchallenge.trajectory.propagator import propagate_cowell, dsat_dt
from solarchallenge.trajectory.perturbations import constant_accel_wrapper

J2000 = datetime(2000, 1, 1, 11, 58, 55, 816) #J2000 in UTC format
AU = 149_597_870.7 #Astronomical Unit in km

#https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/17_frames_and_coordinate_systems.pdf
#https://nssdc.gsfc.nasa.gov/planetary/factsheet/fact_notes.html


"""
### Test for iss
from poliastro.examples import iss
from poliastro.twobody.propagation import CowellPropagator
from poliastro.twobody.sampling import EpochsArray
from astropy import units as u

accel = 2E-2

def constant_accel_factory(accel):
    def constant_accel(t0, u, k):
        v = u[3:]
        norm_v = (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) ** 0.5
        return accel * v / norm_v

    return constant_accel

def f(t0, state, k):
    du_kep = dsat_dt(t0, state, k)
    ax, ay, az = constant_accel_factory(accel)(t0, state, k)
    du_ad = np.array([0, 0, 0, ax, ay, az])

    return du_kep + du_ad

def state_to_vector(ss):
    r, v = ss.rv()
    x, y, z = r.to_value(u.km)
    vx, vy, vz = v.to_value(u.km / u.s)
    return np.array([x, y, z, vx, vy, vz])

r = [  859.07256, -4137.20368,  5295.56871]
v = [7.37289205, 2.08223573, 0.43999979]

iss_orbit = Orbit.from_vector(r,v,J2000, Earth)
print(vars(iss._state))

rtol = 1e-13
full_periods = 2
tf = 100 * u.s#(2 * full_periods + 1) * iss.period / 2 * 0
u0 = state_to_vector(iss)

iss_f_num = iss.propagate(tf, method=CowellPropagator(rtol=rtol, f=f))
print(iss_f_num.r, iss_f_num.v, tf)
print(propagate_cowell(iss_orbit, tofs=100, perturbations=[constant_accel_wrapper(accel)]))

a = constant_accel_wrapper
print(a)
"""

from datetime import timedelta
from src.solarchallenge.model import Model

#earth_orbit = Orbit.from_coe(a=1.00000011*AU, ecc=0.01671022, inc = 0.00005, raan = -11.26064, argp=102.94719, nu = 356.907,
#                             epoch=J2000, attractor=Sun)

earth_orbit = Orbit.from_coe(a=1.00000011*AU, ecc=0.0, inc = 0.0, raan = 0, argp=0, nu = 0,
                             epoch=J2000, attractor=Sun)


#r = [859.07256, -4137.20368,  5295.56871]
#v = [7.37289205, 2.08223573, 0.43999979]
#iss_orbit = Orbit.from_vector(r,v,J2000, Earth)

iss_orbit = Orbit.from_coe(a = 6771, ecc = 0.0, inc=60, raan=0, argp=0, nu = 180, epoch=J2000, attractor=Earth)

engine = Model()
engine.set_orbit_body(earth_orbit)
engine.set_orbit_solar_panel(iss_orbit)
engine.propagate_orbit(start = J2000, end = J2000 + timedelta(0.2), perturbations=[constant_accel_wrapper(-0E-5)])

###################
# Solar Panel stuff

from solarchallenge.solarpanel.solarpanel import SolarPanel

solar_panel = {
    'a': 1,
    'n_s': 1,
    'i_sc': 0.506,
    'v_oc': 2.667,
    'i_mp': 0.487,
    'v_mp': 2.371,
    'alpha_i_sc': 0.32/1000.0,
    'alpha_v_oc': -6.0/1000.0,
    'alpha_v_mp': -6.1/1000.0,
    'alpha_i_mp': 0.28/1000.0,
    'tref': 28,
    'solar_irr_ref': 1367,
    'eff': 1
}

panel = SolarPanel.from_dict(solar_panel)
engine.set_solar_panel(panel)

solar_add = {"area": 1,
             "mass": 1,
             "cp": 1000,
             "emittance_back": 1.0,
             "emittance_front": 1.0,
             "absorptance_back": 1.0,
             "absorptance_front": 1.0,
}

engine.compute_power(thermal_model=True, a=solar_add)

TCELL = 20
RADIANCE = 1300

#print(panel.compute_angle_solar_incidence(r_earth=[300000000,0,0], r_panel=[30.0,0,0]))
#print(panel.compute_power(T = 28, angle_incidence=75, V = 2))

#sol = compute_IV_PindadoCubas(G=RADIANCE, T=TCELL, solar_panel=m, method='constant', V = 2.6)

# Sun distance
# Temperature

# (Initial Tcell, method, rtol, atol, params for albedo, etc...)
