import meep as mp
from matplotlib import pyplot as plt
import numpy as np
from IPython.display import Video
#%matplotlib inline # IPython command used in Jupyter
from pathlib import Path

savedir = Path(__file__).parent / 'results_2'


# Define simulation domain and grid resolution
cell = mp.Vector3(16,16,0)
resolution = 10

# Define waveguide
geometry = [mp.Block(mp.Vector3(12,1,mp.inf),
                center=mp.Vector3(-2.5,-3.5),
                material=mp.Medium(epsilon=12)),
            mp.Block(mp.Vector3(1,12,mp.inf),
                center=mp.Vector3(3.5,2),
                material=mp.Medium(epsilon=12))]
pml_layers = [mp.PML(1.0)]

# Assign boundary
pml_layers = [mp.PML(1.0)]

# Define source
sources = [mp.Source(mp.ContinuousSource(wavelength=2*(11**0.5), width=20),
            component=mp.Ez,
            center=mp.Vector3(-7,-3.5),
            size=mp.Vector3(0,1))]

# Define simulation
sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)

# Visualise
f_geom1 = plt.figure(dpi=150)
sim.plot2D(ax = f_geom1.gca())

# Animate
f_anim = plt.figure(dpi=150)
Animate = mp.Animate2D(sim, fields=mp.Ez, f=f_anim, realtime=False, normalize=True)
sim.run(mp.at_every(0.5,Animate),until=100)

# Save animation
filename = savedir / "bent_waveguide_poor.mp4"
fps = 10
Animate.to_mp4(fps,filename)
Video(filename)

# Redo with better design
sim.reset_meep()
cell = mp.Vector3(16,40,0)
geometry = [mp.Block(mp.Vector3(12,1,mp.inf),
                center=mp.Vector3(-2.5,-3.5),
                material=mp.Medium(epsilon=12)),
            mp.Block(mp.Vector3(1,42,mp.inf),
                center=mp.Vector3(3,5,17),
                material=mp.Medium(epsilon=12))]
# Update cell
sim.cell_size = cell
# Update geometry and shift coordinate system
sim.geometry = geometry
sim.geometry_center = mp.Vector3(0,12,0)

# Re-run until time = 200
sim.run(until=200)

# Plot new design
f_geom2 = plt.figure(dpi=150)
sim.plot2D(ax = f_geom2.gca(), fields=mp.Ez)
# Save figure
f_geom2.savefig(savedir / 'bent_waveguide.png', format='png')

# Generate a x-t slice plot instead
vals = []
def get_slice(sim):
    # get_array can take out a part of the cell
    vals.append(sim.get_array(center=mp.Vector3(0,-3.5),
        size=mp.Vector3(16,0), component=mp.Ez))
sim.reset_meep()
sim.run(mp.at_beginning(mp.output_epsilon),
        mp.at_every(0.6,get_slice),
        until=200)

# x-y-plot
f_slice = plt.figure(dpi=150)
plt.imshow(vals, interpolation='spline36', cmap='RdBu')
plt.axis('off')

# Finish and show
print("Done! saved animation to "+savedir.as_posix())
plt.show()