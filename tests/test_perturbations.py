import numpy as np
from solarchallenge.trajectory.orbit import Orbit, rv2coe
from solarchallenge.bodies.bodies import Earth, EARTH
from solarchallenge.propagator.propagator import propagate_cowell
from solarchallenge.propagator.perturbations import J2_perturbation_wrapper, J3_perturbation_wrapper
from datetime import datetime, timedelta

# Perturbation tests are based on the tests conducted with poliastro
# https://github.com/poliastro/poliastro/blob/main/tests/tests_twobody/test_perturbations.py

def test_J2_propagation_Earth():

    r0 = np.array([-2384.46, 5729.01, 3050.46])  # km
    v0 = np.array([-7.36138, -2.98997, 1.64354])  # km/s

    orbit = Orbit.from_vector(r = r0, v =  v0, epoch=None, attractor=Earth)
    __, __, __, raan0, argp0, __ = rv2coe(r=r0, v=v0, mu=Earth.mu)

    #TODO
    tf = 48.0*3600.0
    J2 = 0.00108263
    R = 6378.1366

    r, v = propagate_cowell(orbit=orbit, tf=tf, perturbations=[J2_perturbation_wrapper(J2, R)])
    __, __, __, raan, argp, __ = rv2coe(r=r, v=v, mu=Earth.mu)

    raan_variation_rate = (raan - raan0) / tf * 3600 # [deg/h]
    argp_variation_rate = (argp - argp0) / tf * 3600 # [deg/h]

    assert np.isclose(raan_variation_rate, -0.172, rtol = 1E-2)
    assert np.isclose(argp_variation_rate, 0.282, rtol = 1E-2)


def test_J3_propagation_Earth():

    a0 = 8970.667
    ecc0 = 0.25
    raan0 = 1.047 / np.pi * 180.0
    nu0 = 0.0
    argp0 = 1.0 / np.pi * 180.0
    inc0 = 0.2618 / np.pi * 180.0

    #From PoliAstro
    r_J2_expected = [ 3882.05385499, -9440.99634518, -2727.32447055]
    r_J3_expected = [ 3865.22644606, -9447.01398762, -2726.87798957]

    orbit = Orbit.from_coe(a=a0, ecc=ecc0, inc=inc0, raan=raan0, argp=argp0, nu=nu0, epoch=None, attractor=Earth)

    #TODO
    tf = 10.0 * 86400.0
    J2 = 0.00108263
    J3 = -2.5326613168e-6
    R = 6378.1366

    t_list = np.linspace(0, tf, 1000)
    r_J2, v_J2 = propagate_cowell(orbit=orbit, tf=tf, t_eval=t_list, perturbations=[J2_perturbation_wrapper(J2, R)])
    r_J3, v_J3 = propagate_cowell(orbit=orbit, tf=tf, t_eval=t_list, perturbations=[J2_perturbation_wrapper(J2, R), J3_perturbation_wrapper(J3, R)])

    assert np.isclose(r_J2[0,-1], r_J2_expected[0], rtol = 1E-2)
    assert np.isclose(r_J2[1,-1], r_J2_expected[1], rtol = 1E-3)
    assert np.isclose(r_J2[2,-1], r_J2_expected[2], rtol = 1E-3)

    assert np.isclose(r_J3[0,-1], r_J3_expected[0], rtol = 1E-2)
    assert np.isclose(r_J3[1,-1], r_J3_expected[1], rtol = 1E-3)
    assert np.isclose(r_J3[2,-1], r_J3_expected[2], rtol = 1E-3)