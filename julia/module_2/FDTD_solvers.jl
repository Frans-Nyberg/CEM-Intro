module FDTDSolvers

# Solves the wave equation using central difference
# Er0 and Er1 are initial conditions at time t=0 and t=dt
# C is Courant number
# Returns a histogram of solutions at each time step
function explicit_1D_wave(Er0::Vector, Er1::Vector, M, N, C2) :: Matrix
    Er = hcat(Er0, Er1)
    Ern = zeros(M,1)
    Epp = copy(Er0);    Dir_bound(Epp)
    Ep = copy(Er1);     Dir_bound(Ep)
    for n = 1:N
        # Step
        for m = 2:M-1
            Ern[m] = 2.0*Ep[m] - Epp[m] + C2*(
                Ep[m+1] - 2.0*Ep[m] + Ep[m-1]
            )
        end
        # Update
        Epp = copy(Ep)
        Ep = copy(Ern)
        # Store
        Er = hcat(Er, Ern)
    end
    return Er
end

function Dir_bound(E::Vector)
    E[1] = 0.
    E[end] = 0.
end

gaussian(x,sig,xpeak) = exp.( -(x.-xpeak).^2 ./ (2*sig^2) )

end