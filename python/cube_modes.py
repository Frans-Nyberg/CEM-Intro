from module_2.modes import Modes
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import meep as mp
from random import random, seed
from mayavi import mlab

savedir = Path(__file__).parent / "results_2"

## Tests and diagnosing: a 2D problem
def test_harmonics(Init_val):
    # Initialise
    N = int(1e3)
    lres = 1e3
    dims = [1,1,0]
    cavity_2d = Modes(N, lres)
    field = mp.Ey
    cavity_2d.init_sim(dims, Init_val)
    anim_2d = mp.Animate2D(cavity_2d.sim, fields=field,
                realtime=False, normalize=True, labels=True)
    # Sampling
    sample_points = [mp.Vector3(cavity_2d.L/3, 0, 0)]
    samples = []
    def sample(sim):
        samples_t = []
        for pos in sample_points:
            si = sim.get_field_point(field, pos).real
            samples_t.append(si)
        samples.append(samples_t)
    cavity_2d.sim.run(
        mp.at_every(10*cavity_2d.dt,anim_2d),
        mp.after_time(100*cavity_2d.dt,sample),until=cavity_2d.t1)
    # Plots
    K = len(samples)
    tK = np.linspace(0, cavity_2d.t1, K)
    fig_2d, (ax_t2d, ax_f2d) = plt.subplots(2,1)
    fig_2d.subplots_adjust(hspace=.5)
    ax_t2d.plot(tK, samples)
    ax_t2d.set_xlabel("time")
    ax_t2d.set_title("Ey")
    # Spectrum
    f, amplitudes = cavity_2d.spectrum(samples)
    Kp = round(K/200)
    amplitude = amplitudes[0]
    ax_f2d.plot(f[0:Kp], amplitude[0:Kp], label="amplitude")
    ax_f2d.vlines(cavity_2d.lowest_mode[Init_val], min(amplitude), max(amplitude),
                linestyles='dashed', label="fundamental mode")
    ax_f2d.set_xlabel("frequency")
    ax_f2d.legend()
    # Save plots
    animpath_2d = savedir / ("cavity_2d_"+Init_val.name+".mp4")
    fps = 10
    anim_2d.to_mp4(fps,animpath_2d)
    samplepath_2d = savedir / ("sample_2d_"+Init_val.name+".png")
    fig_2d.savefig(samplepath_2d, format="png")

def find_cube_modes():
    # Initialise
    N = int(5e3)
    lres = 1e2
    dims = [1,1,1]
    cavity_3d = Modes(N, lres)
    cavity_3d.init_sim(dims, Modes.Init_vals.Random)
    # Run and sample
    seed(732)
    rp = lambda: (random() - cavity_3d.L/2)
    sample_points = [mp.Vector3(cavity_3d.L/3**(1/2), cavity_3d.L/3**(1/2), cavity_3d.L/3**(1/2)), 
                    mp.Vector3(rp(), rp(), rp())]
    samples = []
    def sample(sim):
        samples_t = []
#        for comp in [mp.Ex, mp.Ey, mp.Ez]:
        for comp in [mp.Ez]:
            for pos in sample_points:
                si = sim.get_field_point(comp, pos).real
                samples_t.append(si)
            samples.append(samples_t)
    cavity_3d.sim.run(mp.after_time(100*cavity_3d.dt,sample), until=cavity_3d.t1)
    # Initialise plots
    K = len(samples[0])
    tK = np.linspace(0, cavity_3d.t1, K)
    fig_2d, (ax_t2d, ax_f2d) = plt.subplots(2,1)
    fig_2d.subplots_adjust(hspace=.5)
    # Plot field
    for sample in samples:
        ax_t2d.plot(tK, samples)
    ax_t2d.set_xlabel("time")
    ax_t2d.set_title("Ez")
    # Save plots
    samplepath_2d = savedir / ("sample_2d_"+Modes.Init_vals.Random.name+".png")
    fig_2d.savefig(samplepath_2d, format="png")

#mlab.contour3d(Ex, contours=10)
#mlab.colorbar()

# Reference
K = 3
#v = an_modes(K,K,K)
#ax_spc.plot(0,0)
#ax_spc.vlines(v, c_[0], c_[-1])

if __name__ == "__main__":
#    test_harmonics(Modes.Init_vals.Fundamental)
#    test_harmonics(Modes.Init_vals.FirstHarmonic)
    test_harmonics(Modes.Init_vals.Random)
#    find_cube_modes()
    print("done!")