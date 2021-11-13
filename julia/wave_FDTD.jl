include("module_2/FDTD_solvers.jl")
import .FDTDSolvers: explicit_1D_wave, gaussian
using Plots

savedir = joinpath(pwd(), "results_2")

# domain size
a = 1
# standard deviation
sig = 0.05
# grid points
M = 501

c = 1
dx = a/(M-1)
xm = LinRange(0, a, M)
dt = 0.02
N = 2
Cour = (c*dt/dx)^2

# initialise
Er0 = gaussian(xm, sig, a/2)
Er1 = gaussian(xm, sig, a/2+c*dt)

## Solve
Er = explicit_1D_wave(Er0 ,Er1, M, N, Cour)

## Plot
plot(xm, Er)