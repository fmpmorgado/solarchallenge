import numpy as np
import numpy.typing as npt
from solarchallenge.trajectory.propagator import propagate_cowell
from solarchallenge.solarpanel.solarpanel import SolarPanel
from solarchallenge.trajectory.orbit import Orbit
from solarchallenge.solarpanel.thermal import energy_solver
from solarchallenge.solarpanel.utils import compute_solar_radiance, compute_angle_solar_incidence
from dataclasses import dataclass

@dataclass
class OrbitData:
    """
    Class to store the Orbit data after propagation (inherited by the Model class)

    Parameters:
    ----------
    time : np.array[float]
        Correspondent time of stored position and velocity in UTC format
    r_panel : np.array[float]
        Solar panel position array in x, y and z direction - J2000 frame [km]
    v_panel : np.array[float]
        Solar panel velocity array in x, y and z direction - J2000 frame [km/s]
    r_planet : np.array[float]
        Earth position array in x, y and z direction - EME2000 frame [km]
    v_planet : np.array[float]
        Earth velocity array in x, y and z direction - EME2000 frame [km/s]    
    """

    time: npt.ArrayLike
    r_panel: npt.ArrayLike
    v_panel: npt.ArrayLike
    r_planet: npt.ArrayLike
    v_planet: npt.ArrayLike

    def store_orbit_data(self, time: npt.ArrayLike, r_panel: npt.ArrayLike,
                         v_panel: npt.ArrayLike, r_planet:npt.ArrayLike, v_planet: npt.ArrayLike) -> None:
        """
        Stores the orbit data in the respective object variables and checks for consistency.

        Parameters:
        ----------
        time : np.array[float]
            Correspondent time of stored position and velocity in UTC format
        r_panel : np.array[float]
            Solar panel position array in x, y and z direction - J2000 frame [km]
        v_panel : np.array[float]
            Solar panel velocity array in x, y and z direction - J2000 frame [km/s]
        r_planet : np.array[float]
            Earth position array in x, y and z direction - EME2000 frame [km]
        v_planet : np.array[float]
            Earth velocity array in x, y and z direction - EME2000 frame [km/s] 
        """

        # Check if lists are same size
        if len(time) != len(r_panel) or \
           len(time) != len(v_panel) or \
           len(time) != len(r_planet) or \
           len(time) != len(v_planet):
            raise ValueError("Lists need to have the same size.")

        # Store the orbit data
        self.time = time
        self.r_panel = r_panel
        self.v_panel = v_panel
        self.r_planet = r_planet
        self.v_planet = v_planet


@dataclass
class HeatData:
    q_earth_ir: npt.ArrayLike
    q_earth_alb: npt.ArrayLike
    q_solar_front: npt.ArrayLike
    q_solar_back: npt.ArrayLike
    q_sat_rad: npt.ArrayLike

    def store_heat_data(self, q_earth_ir, q_earth_alb, q_solar_front, q_solar_back, q_sat_rad):
        self.q_earth_ir = q_earth_ir
        self.q_earth_alb = q_earth_alb
        self.q_solar_front = q_solar_front
        self.q_solar_back = q_solar_back
        self.q_sat_rad = q_sat_rad
        

@dataclass
class PowerData:
    voltage: npt.ArrayLike
    current: npt.ArrayLike
    power: npt.ArrayLike
    temperature: npt.ArrayLike

    def store_power_data(self, voltage, current, power, temperature):
        self.voltage = voltage
        self.current = current
        self.power = power
        self.temperature = temperature


class Model(OrbitData, PowerData, HeatData):
    def __init__ (self):
        self.orbit_panel: Orbit
        self.orbit_body: Orbit
        self.panel: SolarPanel

    def set_orbit_solar_panel(self, orbit_panel):
        self.orbit_panel = orbit_panel

    def set_orbit_body(self, orbit_body):
        self.orbit_body = orbit_body

    def set_solar_panel(self, panel):
        self.panel = panel

    def propagate_orbit(self, start, end, npoints = 10000, method = 'cowell', ivp_method = 'DOP853', rtol = 1E-11, atol = 1E-12, perturbations = None):

        #Verify if orbits are stored:
        if not self.orbit_panel:
            raise Exception("Need to set the orbit of the solar panel.")
        if not self.orbit_body:
            raise Exception("Need to set the orbit of the planet (ex: Earth orbit).")

        #Verify if end date is after start date:
        if end < start:
            raise ValueError("Start date needs to be before end date")
        
        #Verify if start date is not before orbit epochs
        if start < self.orbit_panel.epoch or start < self.orbit_body.epoch:
            raise ValueError("Start date needs to be after orbit epoch")
 
        t0_body = (start - self.orbit_body.epoch).total_seconds()
        tf_body = (end - self.orbit_body.epoch).total_seconds()

        t0_panel = (start - self.orbit_panel.epoch).total_seconds()
        tf_panel = (end - self.orbit_panel.epoch).total_seconds()

        t_list_body = np.linspace(t0_body, tf_body, npoints)
        t_list_panel = np.linspace(t0_panel, tf_panel, npoints)

        r_planet, v_planet = propagate_cowell(orbit = self.orbit_body, tf = tf_body, t_eval = t_list_body)
        r_panel, v_panel = propagate_cowell(orbit = self.orbit_panel, tf = tf_panel,
                                            t_eval = t_list_panel, perturbations=perturbations)

        t_list_UTC = np.linspace(start.timestamp(), end.timestamp(), npoints)

        r_planet = np.swapaxes(r_planet,0,1)
        v_planet = np.swapaxes(v_planet,0,1)
        r_panel = np.swapaxes(r_panel,0,1)
        v_panel = np.swapaxes(v_panel,0,1)

        self.store_orbit_data(time=t_list_UTC, r_panel=r_panel, v_panel=v_panel, r_planet=r_planet, v_planet=v_planet) 

    def compute_power(self, temperature=28, voltage=2.6, method_power='pindado_cubas',
                      maximum_power=False, thermal_model = False, a=None):

        if not self.panel:
            raise Exception("Need to set the solar panel.")

        if thermal_model == True and not a:
            raise Exception("XXXXX parameters need to be specified to perform thermal modelling.")

        if maximum_power is True:
            voltage = None

        #Compute angle of incidence
        theta = compute_angle_solar_incidence(self.r_planet, self.r_panel)
        
        #Compute solar irradiance
        solar_irr = compute_solar_radiance(self.r_planet, self.r_panel)

        if thermal_model is False:
            power, cur, vol = self.panel.compute_power(method = method_power, solar_irr=solar_irr,
                                                angle_incidence = theta, temperature=temperature, voltage=voltage)

            self.store_power_data(voltage=vol, current=cur, power=power, temperature=np.ones(power.shape) * temperature)
    
        else: #Thermal model is active
            power, cur, vol, temp, heat = energy_solver(panel=self.panel, time=self.time,
                                                                       r_panel=self.r_panel,
                                                                       r_planet=self.r_planet,
                                                                       solar_irr=solar_irr, theta=theta,
                                                                       heat_parameters=a, voltage=voltage)
    
            self.store_power_data(voltage=vol, current=cur, power=power, temperature=temp)
            self.store_heat_data(**heat)
        """
        import plotly.express as px
        import pandas as pd

        df = pd.DataFrame(dict(
        x = self.time,
        y = T
        ))

        df = pd.DataFrame(Q)
        df['x'] = self.time

        fig = px.line(df, x = 'x', y = df.columns[:-1], title='Life expectancy in Canada')
        fig.show()

        """

"""
import plotly.graph_objs as go

def ms(x, y, z, radius, resolution=20):
    u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
    X = radius * np.cos(u)*np.sin(v) + x
    Y = radius * np.sin(u)*np.sin(v) + y
    Z = radius * np.cos(v) + z
    return (X, Y, Z)

def plot_orbit(r):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter3d(
            x=r[:,0],
            y=r[:,1],
            z=r[:,2],
            mode='lines',
        ))
    
    (X, Y, Z) = ms(0,0,0,6000)
    fig.add_trace(
        go.Surface(
            x=X,
            y=Y,
            z=Z,
        ))


    #fig.show(renderer="png")
#from plotly.offline import iplot, init_notebook_mode
#init_notebook_mode()
#iplot(fig)
"""