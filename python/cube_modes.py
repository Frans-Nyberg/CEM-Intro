from module_2.modes import Modes
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import meep as mp
from random import random, seed

savedir = Path(__file__).parent / "results_2"

## The 3D cubic box
def find_cube_modes():
    # Initialise
    N = int(1e4)
    lres = 3e1
    dims = [1,1,1]
    Courant = 1/3**(1/2)
    cavity_3d = Modes(N, lres, Courant)
    cavity_3d.init_sim(dims, Modes.Init_vals.Random)
    # Run and sample
    seed()
    rp = lambda: (random() - cavity_3d.L/2)
    sample_points = [mp.Vector3(cavity_3d.L*3**(1/2)/10, cavity_3d.L*3**(1/2)/10, cavity_3d.L*3**(1/2)/10), 
                    mp.Vector3(rp(), rp(), rp())]
    point_to_str = lambda i: "(%.2f,%.2f,%.2f)" %(sample_points[i].x, sample_points[i].y, sample_points[i].z)
    sample_names = ["Manual pick: "+point_to_str(0), "Random pick: "+point_to_str(1)]
    sample_linestyle = ["-g", "--r"]
    samples = []
    def sample(sim):
        samples_t = []
        for comp in [mp.Ez]:    # symmetry, only need one component
            for pos in sample_points:
                si = sim.get_field_point(comp, pos).real
                samples_t.append(si)
            samples.append(samples_t)
    cavity_3d.sim.run(mp.after_time(cavity_3d.dt,sample), until=cavity_3d.t1)
    # Initialise plots
    K = len(samples)
    tK = np.linspace(0, cavity_3d.t1, K)
    fig_3d, (ax_t3d, ax_f3d) = plt.subplots(2,1)
    fig_3d.subplots_adjust(hspace=.5)
    # Plot field
    Kp = round(K/50)
    for i, sample in enumerate(np.transpose(samples)):
        ax_t3d.plot(tK[0:Kp], sample[0:Kp], sample_linestyle[i], label=sample_names[i])
    ax_t3d.set_xlabel("time")
    ax_t3d.set_title("Ez")
    ax_t3d.legend()
    # Spectrum
    f, amplitudes = cavity_3d.spectrum(samples)
    f_ans = cavity_3d.an_modes(4,5,5)
    for f_an in f_ans:
        # Analytical spectrum
        if f_an <= f[Kp]:
            ax_f3d.vlines(f_an, min(amplitudes[0]), max(amplitudes[0][0:Kp]),  linestyles='dashed')
    for i, amplitude in enumerate(amplitudes):
        # Numerical spectrum
        ax_f3d.plot(f[0:Kp], amplitude[0:Kp], sample_linestyle[i], label=sample_names[i])
    ax_f3d.set_xlabel("frequency")
    ax_f3d.set_title("amplitude")
    ax_f3d.legend()
    # Save plots
    samplepath_2d = savedir / ("sample_3d_"+Modes.Init_vals.Random.name+".png")
    fig_3d.savefig(samplepath_2d, format="png")

## Tests and diagnosing: a 2D problem
def test_harmonics(Init_val, N, lres, sfac=100, Kfac=200):
    # Initialise
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
        mp.at_every(100*cavity_2d.dt,anim_2d),
        mp.after_time(sfac*cavity_2d.dt,sample),until=cavity_2d.t1)
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
    Kp = round(K/Kfac)
    amplitude = amplitudes[0]
    ax_f2d.plot(f[0:Kp], amplitude[0:Kp], label="amplitude")
    ax_f2d.vlines(cavity_2d.lowest_mode[Init_val], min(amplitude), max(amplitude[0:Kp]),
                linestyles='dashed', label="fundamental mode")
    ax_f2d.set_xlabel("frequency")
    ax_f2d.legend()
    # Save plots
    animpath_2d = savedir / ("cavity_2d_"+Init_val.name+".mp4")
    fps = 10
    anim_2d.to_mp4(fps,animpath_2d)
    samplepath_2d = savedir / ("sample_2d_"+Init_val.name+".png")
    fig_2d.savefig(samplepath_2d, format="png")

if __name__ == "__main__":
#    test_harmonics(Modes.Init_vals.Fundamental, int(5e3), 1e2)
#    test_harmonics(Modes.Init_vals.FirstHarmonic, int(5e3), 1e2)
#    test_harmonics(Modes.Init_vals.Random, int(1e4), 3e1, 1, 50)
    find_cube_modes()
    print("done!")