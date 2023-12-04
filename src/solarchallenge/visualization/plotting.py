"""
Module to plot orbits and computed data

Caution: This module is for demo purposes
"""
from solarchallenge.trajectory.orbit import rv2coe
from solarchallenge.solarpanel.solarpanel import SolarPanel

import plotly.graph_objs as go
from plotly.subplots import make_subplots
import numpy as np
from datetime import datetime

def plot_orbit(r, radius_body):
    """
    Plot orbit around attractor
    """
    
    pos_earth = np.array([0.0, 0.0, 0.0])

    def sphere(x, y, z, radius, resolution=20):
        u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
        X = radius * np.cos(u)*np.sin(v) + x
        Y = radius * np.sin(u)*np.sin(v) + y
        Z = radius * np.cos(v) + z
        return (X, Y, Z)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter3d(
            x=r[:,0],
            y=r[:,1],
            z=r[:,2],
            mode='lines'
        ))
    
    (X, Y, Z) = sphere(pos_earth[0], pos_earth[1], pos_earth[2], radius_body)
    
    fig.add_trace(
        go.Surface(
            x=X,
            y=Y,
            z=Z,
            showscale=False
        ))  
    
    fig.show()


def plot_coe(r, v, time, mu):
    """
    Plot classical orbital elements
    """
    fig = make_subplots(rows=2, cols=3,
                        subplot_titles=("Semi-major axis",
                                        "Eccentricity",
                                        "Inclination",
                                        "RAAN",
                                        "Argument of periapsis",
                                        "True anomaly"))
    #Timestamp to datetime
    time = time.astype(int)
    time = np.array([datetime.fromtimestamp(t) for t in time])
    
    #Compute orbit classical orbital elements
    a = np.zeros(time.shape)
    ecc = np.zeros(time.shape)
    inc = np.zeros(time.shape)
    raan = np.zeros(time.shape)
    argp = np.zeros(time.shape)
    nu = np.zeros(time.shape)

    for i, (r_i, v_i) in enumerate(zip(r, v)):
        a[i], ecc[i], inc[i], raan[i], argp[i], nu[i] = rv2coe(r_i, v_i, mu)

    #Rounding the values
    a = np.round(a, 5)
    ecc = np.round(ecc, 8)
    inc = np.round(inc, 8)
    raan = np.round(raan, 8)
    argp = np.round(argp, 8)
    nu = np.round(nu, 8)

    fig.add_trace(go.Scatter(x=time, y=a), row=1, col=1)
    fig.add_trace(go.Scatter(x=time, y=ecc), row=1, col=2)
    fig.add_trace(go.Scatter(x=time, y=inc), row=1, col=3)
    fig.add_trace(go.Scatter(x=time, y=raan), row=2, col=1)
    fig.add_trace(go.Scatter(x=time, y=argp), row=2, col=2)
    fig.add_trace(go.Scatter(x=time, y=nu), row=2, col=3)

    fig.update_yaxes(title_text="km",  row=1, col=1)
    fig.update_yaxes(title_text="deg", row=1, col=3)
    fig.update_yaxes(title_text="deg", row=2, col=1)
    fig.update_yaxes(title_text="deg", row=2, col=2)
    fig.update_yaxes(title_text="deg", row=2, col=3)

    fig.update_xaxes(title_text="Date", row=2, col=1)
    fig.update_xaxes(title_text="Date", row=2, col=2)
    fig.update_xaxes(title_text="Date", row=2, col=3)
    
    fig.show()

def plot_IV(panel: SolarPanel, solar_irr, temperature):
    """
    Plot IV of solar panel
    """

    #Check v_mp
    v_mp = panel.v_mp
    v_max = 1.5 * v_mp

    v_array = np.linspace(0, v_max, 200)

    fig = make_subplots(rows=1, cols=2, subplot_titles=("Voltage-Current", "Voltage-Power"))

    try:
        for t in temperature:
            p, i, __ = panel.compute_power(solar_irr=solar_irr, temperature=t*np.ones(v_array.shape), voltage=v_array)

            fig.add_trace(go.Scatter(x=v_array, y=i,
                    mode='lines',
                    name=f"Temperature: {t} ºC"),
                    row=1, col=1)
            

            fig.add_trace(go.Scatter(x=v_array, y=p,
                    mode='lines',
                    name=f"Temperature: {t} ºC"),
                    row=1, col=2) 

    except:
        pass

    try:
        for g in solar_irr:

            p, i, __ = panel.compute_power(solar_irr=g * np.ones(v_array.shape), temperature=temperature, voltage=v_array)

            fig.add_trace(go.Scatter(x=v_array, y=i,
                    mode='lines',
                    name=f"Effective solar irradiance: {g} W"),
                    row=1, col=1)
            

            fig.add_trace(go.Scatter(x=v_array, y=p,
                    mode='lines',
                    name=f"Effective solar irradiance: {g} W"),
                    row=1, col=2) 

    except:
        pass


    fig.update_yaxes(title_text="Current [A]", row=1, col=1)
    fig.update_yaxes(title_text="Power [W]", row=1, col=2)

    fig.update_xaxes(title_text="Voltage [V]", row=1, col=1)
    fig.update_xaxes(title_text="Voltage [V]", row=1, col=2)

    fig.show()