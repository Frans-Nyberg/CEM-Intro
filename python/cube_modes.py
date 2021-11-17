from module_2.eigenfuncs import iter_trips
import meep as mp
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from random import random, seed
from mayavi import mlab

savedir = Path(__file__).parent / "results_2"

# Cube length
L = 1
# Light speed in cavity
c = 1
# Time steps
N = 1e3

# Analytical modes
an_f = lambda v: 0.5*c/L * (v[0]**2 + v[1]**2 + v[2]**2)**(1/2)
an_fmin = an_f([1,0,0])
an_modes = lambda N,M,P: \
    [an_f(v) for v in iter_trips(N,M,P)]

## MEEP
# Define simulation domain
cell = mp.Vector3(L,0,0)
# define resolution in terms of number of largest wavelengths
lam_res = 10
res = lam_res * round(c/an_fmin)
# Courant number
C = 0.5

# Initialise simulation
sim = mp.Simulation(cell, res, Courant=C)
# simulation time
t1 = N*C/c/res

# Define random initial values
seed(2370)
make_rand = lambda pos: random() - 0.5
components = [mp.Ex, mp.Ey, mp.Ez]
for comp in components: sim.initialize_field(comp, make_rand)

# Run simulation and get fields
sim.run(until=t1)

# Get fields
Ex = sim.get_efield_x()

# Check on the fields
f_ch, ax_ch = plt.subplots()
sim.plot2D(ax=ax_ch, fields=mp.Ex)
#mlab.contour3d(Ex, contours=10)
#mlab.colorbar()

f_ = np.linspace(0,1)
c_ = np.linspace(0,1)

## Spectral Plot
fig_spc, ax_spc = plt.subplots()

# Reference
K = 3
v = an_modes(K,K,K)
ax_spc.plot(0,0)
ax_spc.vlines(v, c_[0], c_[-1])

if __name__ == "__main__":
#    mlab.show()
    plt.show()
    print("done!")