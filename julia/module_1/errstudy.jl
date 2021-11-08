module ErrStudy
using BenchmarkTools

function find_rounderr(op, maker, op_round, norms, S::Int,
Ms::AbstractArray
) :: Matrix{Float64}
    errM = fill([], length(Ms))
    for (i, M) in enumerate(Ms)
        err_norms = zeros(length(norms))
        for _ in 1:S
            d1 = maker(M);
            d2 = maker(M);
            outp = op(d1, d2)
            d1_round = op_round(d1);
            d2_round = op_round(d2);
            outp_round = op(d1_round, d2_round)
            for (j, norm) in enumerate(norms)
                err_norms[j] += norm(outp_round-outp)
            end
        end
        errM[i] = err_norms/S
    end
    # transpose and convert to matrix
    return hcat(errM...)
end

function time_op(op, maker, Ms)
    tM = zeros(length(Ms))
    for (i, M) in enumerate(Ms)
        d1 = maker(M); d2 = maker(M)
        # $ because local variables
        # 1e9 because nanoseconds
        b = @benchmark $op($1, $d2)
        tM[i] = median(b.times) / 1e9
    end
    return tM
end

end