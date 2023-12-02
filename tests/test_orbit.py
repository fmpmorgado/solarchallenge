import numpy as np
from solarchallenge.trajectory.orbit import Orbit
from solarchallenge.trajectory.orbit import coe2rv, rv2coe
from solarchallenge.bodies.bodies import Earth

def test_convert_coe_to_rv():

    earth_mu = Earth.mu
    p = 11_067.790
    ecc = 0.83285
    inc = 87.87
    raan = 227.89
    argp = 53.38
    nu = 92.335

    expected_r = [6_525.344, 6_861.535, 6_449.125]
    expected_v = [4.902276, 5.533124, -1.975709]

    r, v = coe2rv(a=p / (1 - ecc**2),
                  ecc=ecc,
                  inc=inc,
                  raan=raan,
                  argp=argp,
                  nu=nu,
                  mu=earth_mu
                  )

    assert np.isclose(r, expected_r, rtol=1e-5).all()
    assert np.isclose(v, expected_v, rtol=1e-5).all()


def test_convert_rv_to_coe():

    earth_mu = Earth.mu
    p_expected = 11_067.790
    ecc_expected = 0.83285
    inc_expected = 87.87
    raan_expected = 227.89
    argp_expected = 53.38
    nu_expected = 92.335

    r = [6_525.344, 6_861.535, 6_449.125]
    v = [4.902276, 5.533124, -1.975709]

    a, ecc, inc, raan, argp, nu = rv2coe(r, v, mu=earth_mu)

    assert np.isclose(a * (1 - ecc**2), p_expected, rtol=1e-5)
    assert np.isclose(ecc, ecc_expected, rtol=1e-5)
    assert np.isclose(inc, inc_expected, rtol=1e-5)
    assert np.isclose(raan, raan_expected, rtol=1e-5)
    assert np.isclose(argp, argp_expected, rtol=1e-5)
    assert np.isclose(nu, nu_expected, rtol=1e-5)


def test_rv_initialization_equal_to_coe():
    
    #classical orbital elements
    p = 11_067.790
    ecc = 0.83285
    inc = 87.87
    raan = 227.89
    argp = 53.38
    nu = 92.335

    #rv elements
    r = [6_525.344, 6_861.535, 6_449.125]
    v = [4.902276, 5.533124, -1.975709]

    coe_orbit = Orbit.from_coe(a=p / (1 - ecc**2),
                                ecc=ecc,
                                inc=inc,
                                raan=raan,
                                argp=argp,
                                nu=nu,
                                epoch=0,
                                attractor=Earth)
    
    rv_orbit = Orbit.from_vector(r=r, v=v, epoch=0, attractor=Earth)

    assert np.isclose(coe_orbit.r, rv_orbit.r, rtol=1e-5).all()
    assert np.isclose(coe_orbit.v, rv_orbit.v, rtol=1e-5).all()
