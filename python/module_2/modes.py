import meep as mp
from random import random, seed
from enum import Enum
from math import cos, sin, pi
import numpy as np

# Class that helps making square or cubic cavity simulations
class Modes:
    # Cube length
    L = 1
    # Light speed in cavity
    c = 1
    # Initial values
    Init_vals = Enum('Initial_values', 'Fundamental FirstHarmonic Random')
    # Lowest mode corresponding
    # Initial field maximum
    E0 = 1
    # Simulation object
    sim = None
    # Cells
    cell = None

    def __init__(self, N, lam_res, Courant=0.5):
        # Time steps
        self.N = N
        # Courant number
        self.Cour = Courant
        # define resolution in terms of number of largest wavelengths
        self.res = lam_res * round(self.c/self.an_f([1,0,0]))
        # simulation time
        self.dt = self.Cour/self.c/self.res
        self.t1 = N*self.dt
        self.lowest_mode = {self.Init_vals.Fundamental: self.an_f([1,0,0]),
                            self.Init_vals.FirstHarmonic: self.an_f([2,0,0]),
                            self.Init_vals.Random: self.an_f([1,0,0])
                            }

    # Initialise a simulation problem
    # dimensions dims is [xbool, ybool, zbool]
    def init_sim(self, dims, initial):
        # Define simulation domain
        self.cell = mp.Vector3(self.L*dims[0], self.L*dims[1], self.L*dims[2]) 
        # Initialise simulation
        self.sim = mp.Simulation(self.cell, self.res, Courant=self.Cour)
        components = [mp.Ex, mp.Ey, mp.Ez]
        if initial == self.Init_vals.Random:
            # Define random initial values
            seed(2370)
            def make_rand(pos):
                # Enforce boundary conditions
                wl = lambda comp: abs(comp) < self.L/2
                if wl(pos.x) and wl(pos.y) and wl(pos.z): E = self.E0*(random() - self.L/2)
                else: E = 0
                return E
            # Set random initial fields in simulation
            for comp in components: self.sim.initialize_field(comp, make_rand)        
        if initial == self.Init_vals.Fundamental:
            # Initialise with a fundamental mode (in 2D)
            assert(self.cell.z==0)
            def make_fundamental(pos):
                return self.E0*cos(pi/self.L*pos.x)
            self.sim.initialize_field(mp.Ey, make_fundamental)
        if initial == self.Init_vals.FirstHarmonic:
            # Initialise with a first harmonic mode (in 2D)
            assert(self.cell.z==0)
            def make_fundamental(pos):
                return self.E0*sin(2*pi/self.L*pos.x)
            self.sim.initialize_field(mp.Ey, make_fundamental)

    # Spectral analysis functions
    # samples: [t1:[point 1, point2, ...], t2:[point 1, point2, ...]]
    def spectrum(self, samples):
        K = len(samples)
        hK = round(K/2)-1
        f = np.fft.fftfreq(K, self.dt)
        f_pos = f[0:hK]
#        omega = 2*pi/self.t1 * np.linspace(0, self.N-1, K)
        amplitudes = []
        for samples_pos in np.transpose(samples):
            spectrum = np.fft.fft(samples_pos)
            spectrum_pos = spectrum[0:hK]
            amplitude = np.abs(spectrum_pos)
            amplitudes.append(amplitude)
        return f_pos, amplitudes
        
    # analytical modes
    def an_modes(self, N,M,P):
        return [self.an_f(v) for v in self.iter_trips(N,M,P)]

    def iter_trips(self, M,N,P):
        for m in range(M):
            for n in range(N):
                for p in range(P):
                    if m or n or p:
                        yield (m,n,p)

    def an_f(self, v):
        return 0.5*self.c/self.L * (v[0]**2 + v[1]**2 + v[2]**2)**(1/2)

