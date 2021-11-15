include("module_2/FDTD_solvers.jl")
import .FDTDSolvers: explicit_1D_wave, gaussian
import FFTW: fft
import PyPlot as plt

savedir = joinpath(pwd(), "results_2")

## Simulation to show dispersion

# domain size
a = 1
# ~pulse width
sig = 0.02
# grid points in relation to pulse width
M = round(Int64, 3*a/sig)
# time steps
N = 101

c = 1
dx = a/(M-1)
xm = LinRange(0, a, M)
half(x) = round(Int64,x/2)
dtR(R) = R*dx/c

function propagate(R)::Matrix
    dt = dtR(R)
    # initialise
    xstart = a/5
    Er0 = gaussian(xm, sig, xstart)
    Er1 = gaussian(xm, sig, xstart+c*dt)
    ## Solve
    Er = explicit_1D_wave(Er0 ,Er1, M, N, R^2)
    return Er
end

## Propagate
println("Computing waves")
Er_disp = propagate(0.9)
Er_magic = propagate(1)

## Plot wave
fig_prop, ax = plt.subplots()
function three_plot(Er, linestyle, labl)
    ax.plot(xm, Er[:,1], linestyle, label=labl)
    ax.plot(xm, Er[:,half(N)], linestyle)
    ax.plot(xm, Er[:,end], linestyle)
end
three_plot(Er_magic, "-b", "R = 1")
three_plot(Er_disp, "--r", "R = 0.9")
ax.set_xlabel("z")
ax.set_ylabel("E")
ax.legend()

## Test 1D fft
function test_fft(Er)
    r = Er[:,end]
    s = fft(r)
    s = abs.(s)
    plt.plot(s)
end
fig_test1 = plt.figure()
test_fft(Er_disp)

## Test 2D fft
function test_fft(xm)
    fig, ax = plt.subplots()
    r = sin.(100*pi*xm) * ones(1,length(xm))
    s = fft(r)
    s = abs.(s)
    ons = 1 * (0:M-1)/M
    fx = 0.5/dx * (0:half(M)-1)/half(M)
    c = ax.pcolor(ons, fx, s[1:half(M),:],
        shading="auto")
    return fig
end
fig_test2 = test_fft(xm)

## Test: A non-dispersive spectrum
function test_fft2(E)
    fig, ax = plt.subplots()
    s = fft(E[:,1]) * ones(1,N)
    s = log10.(abs.(s))
    c = ax.pcolor(s[1:half(M),1:half(N)])
    fig.colorbar(c, ax=ax)
    return fig
end
fig_test3 = test_fft2(Er_disp)


## Plot spectra
function plot_spectrum(Er, R, labl)
    fig, ax = plt.subplots()
    # Amplitudes
    spectrum = fft(Er)  # multi-dimensional by default
    logabs = log10.(abs.(spectrum))
    # Axes
    omega = pi/dtR(R) * (0:half(N)-1)/half(N)
    k = pi/dx * (0:half(M)-1)/half(M)
    # Plot
    c = ax.pcolor(omega, k, logabs[1:half(M),1:half(N)],
        vmin=-16, vmax=2, shading="nearest")
    ax.set_xlabel("omega")
    ax.set_ylabel("k")
    ax.set_title(labl)
    fig.colorbar(c, ax=ax)
    return fig
end

fig_disp = plot_spectrum(Er_disp, 0.9, "R = 0.9")
fig_magic = plot_spectrum(Er_magic, 1, "R = 1")

## Save plots
fig_prop.savefig(joinpath(savedir, "gauss_prop.png"), format="png")
fig_disp.savefig(joinpath(savedir, "spctr_disp.png"), format="png")
fig_magic.savefig(joinpath(savedir, "spctr_magic.png"), format="png")

println(string("done! saved in ", savedir))
