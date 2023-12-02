import numpy as np
from datetime import datetime
from src.propagator.propagator import propagate_cowell

class OrbitData:
    def __init__(self):
        self.time = np.array([])
        self.r_panel = np.array([])
        self.v_panel = np.array([])
        self.r_body = np.array([])
        self.v_body = np.array([])
    """
    @property
    def time(self):
        return np.array(self.time)

    @property.setter
    def time(self):


    @property
    def r_panel(self):
        return np.array(self.r_panel)

    @property
    def v_panel(self):
        return np.array(self.v_panel)
    
    @property
    def r_body(self):
        return np.array(self.r_body)

    @property
    def v_body(self):
        return np.array(self.v_body)

    """

    def store_orbit_data(self, time, r_panel, v_panel, r_body, v_body):
        
        # Check if lists are same size
        if len(time) != len(r_panel) or \
           len(time) != len(v_panel) or \
           len(time) != len(r_body) or \
           len(time) != len(v_body):
            raise ValueError("Lists need to have the same size.")

        self.time = time
        self.r_panel = r_panel
        self.v_panel = v_panel
        self.r_body = r_body
        self.v_body = v_body


class Model(OrbitData):
    def __init__ (self):
        self.orbit_panel = None
        self.orbit_body = None
        self.panel = None

    def set_orbit_solar_panel(self, orbit_panel):
        self.orbit_panel = orbit_panel

    def set_orbit_body(self, orbit_body):
        self.orbit_body = orbit_body

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

        r_body, v_body = propagate_cowell(orbit = self.orbit_body, tf = tf_body, t_eval = t_list_body)
        r_panel, v_panel = propagate_cowell(orbit = self.orbit_panel, tf = tf_panel, t_eval = t_list_panel, perturbations=perturbations)

        t_list_UTC = np.linspace(start.timestamp(), end.timestamp(), npoints)

        r_body = np.swapaxes(r_body,0,1)
        v_body = np.swapaxes(v_body,0,1)
        r_panel = np.swapaxes(r_panel,0,1)
        v_panel = np.swapaxes(v_panel,0,1)

        self.store_orbit_data(time=t_list_UTC, r_panel=r_panel, v_panel=v_panel, r_body=r_body, v_body=v_body)
        #plot_orbit(r_panel)

        
import plotly.graph_objs as go

def ms(x, y, z, radius, resolution=20):
    """Return the coordinates for plotting a sphere centered at (x,y,z)"""
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
#    fig.show(renderer="png")
#from plotly.offline import iplot, init_notebook_mode
#init_notebook_mode()
#iplot(fig)