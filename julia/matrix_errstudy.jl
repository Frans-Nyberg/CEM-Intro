include("module_1/errstudy.jl")
import .ErrStudy: find_rounderr, time_op
import LinearAlgebra: norm
using Plots

savedir = joinpath(pwd(), "results_1")

make_mat(M) = rand(Float64, M, M)
matmul(A, B) = A*B

expn = 5:1:11
Ms = 2 .^expn

# Finds round-off error of matrix multiplication
function find_matmul_rounderr() :: Matrix
    mat_round(A) = Array{Float32}(A)
    S = 10
    norms = [
        A -> norm(A, 2),     # Frobenius norm
        A -> norm(A, Inf)   # Max norm
    ]
    return find_rounderr(matmul, make_mat,
        mat_round, norms, S, Ms)
end

# finds time complexity of matrix multiplication
time_matmul() = time_op(matmul, make_mat, Ms)
cmplxn(n, sc) = 2*Ms .^n .*sc/max(Ms...)^n

## Tests

function test_matmul()
    I = [1 0; 0 1]
    A = make_mat(2)
    IA = matmul(I, A)
    @assert(IA[1,2] == A[1,2])
    @assert(IA[2,1] == A[2,1])
end

## Run

println("doing tests")
test_matmul()
println("finding rounding error")
errs = find_matmul_rounderr()
println("finding time complexity (takes a minute)")
times = time_matmul()

## Plot matmul
p1 = plot(Ms, errs[1,:], xaxis=:log,
    yaxis=:log, label="Frobenius-norm",
    legend=:topleft, marker=:circle)
plot!(Ms, errs[2,:], xaxis=:log,
    yaxis=:log, label="max-norm", marker=:x)
plot!(Ms, cmplxn(2, errs[1,end]), linestyle=:dash,
    label="O(n^2)", marker=:star5)
plot!(Ms, cmplxn(1, errs[2,end]), linestyle=:dash,
    label="O(n^1)", marker=:star4)
xaxis!(p1, "rows or columns")
savefig(p1, joinpath(savedir, "rounderr_matmul.png"))

## Plot time
p2 = plot(Ms, times, label="time [s]",
    xaxis=:log, yaxis=:log, legend=:topleft,
    marker=:cross)
plot!(Ms, cmplxn(2, times[end]), label="O(n^2)",
    marker=:star5)
plot!(Ms, cmplxn(3, times[end]), label="O(n^3)",
    marker=:star6)
xaxis!(p2, "rows of columns")
savefig(p2, joinpath(savedir, "time_matmul.png"))

println(string("done! saved in ", savedir))