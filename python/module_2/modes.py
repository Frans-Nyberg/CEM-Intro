from eigenfuncs import iter_trips
import meep as mp
from random import random, seed

# Class that helps making square or cubic cavity simulations
class Modes:
    def __init__(self, N=1e3, Courant=0.5, lam_res=10):
        # Cube length
        self.L = 1
        # Light speed in cavity
        self.c = 1
        # Time steps
        self.N = N
        # Courant number
        self.Cour = Courant
        # define resolution in terms of number of largest wavelengths
        lam_res = 10
        self.res = lam_res * round(self.c/self.an_f([1,0,0]))
        # simulation time
        self.t1 = N*self.Cour/self.c/self.res
        self.sim = None

    def an_f(self, v):
        return 0.5*self.c/self.L * (v[0]**2 + v[1]**2 + v[2]**2)**(1/2)
    
    # analytical modes
    def an_modes(self, N,M,P):
        return [self.an_f(v) for v in iter_trips(N,M,P)]

    # Initialise a simulation problem
    # dimensions dims is [xdim, ydim, zdim]
    def init_sim(self, dims):

        # Define simulation domain
        cell = mp.Vector3(dims[0],dims[1],dims[2])
        # Initialise simulation
        self.sim = mp.Simulation(cell, self.res, Courant=self.Cour)

        # Define random initial values
        seed(2370)
        make_rand = lambda pos: random() - 0.5
        components = [mp.Ex, mp.Ey, mp.Ez]
        for comp in components: self.sim.initialize_field(comp, make_rand)

    # Run the simulation
    def run_sim(self):
        self.sim.run(until=self.t1)