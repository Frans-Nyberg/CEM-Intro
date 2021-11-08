from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from module_1.errstudy import find_rounderr_sweep, time_sweep

savedir = Path(__file__).parent / 'results_1'

# Generates double precision random matrix of values in (0,1) of size M
make_mat = lambda M: np.asmatrix(np.random.rand(M,M), dtype=np.double)
mat_multiply = lambda mat1, mat2: mat1 * mat2

# M sweep range
expn = range(5, 12)
Ms = [2**n for n in expn]

def find_matmul_rounderr():
    # Rounds to single precision
    mat_round = lambda mat: np.asmatrix(mat, dtype=np.single)
    # samples
    S = 10
    norm_arr = [
        lambda mat: np.linalg.norm(mat, 'fro'), # Frobenius norm (euclidean norm for each element)
        lambda mat: np.linalg.norm(mat, np.inf),# matrix Inf norm (max of each abs row sum)
        lambda mat: np.max(np.abs(mat)),        # max norm (just the maximum element)
        lambda mat: np.sum(np.abs(mat)),        # 1-norm (for each element)
        ]
    # transpose to nest the sizes
    return np.transpose(find_rounderr_sweep(mat_multiply, make_mat, mat_round, norm_arr, S, Ms))

def time_matmul():
    S = 100
    return time_sweep(mat_multiply, make_mat, S, Ms)

def plot_matmul_rounderr(err_tup):
    fig, ax = plt.subplots()
    ax.loglog(Ms, err_tup[0], '.-', label="frobenius-norm")
    ax.loglog(Ms, err_tup[1], 'x-', label="inf-norm")
    ax.loglog(Ms, err_tup[2], '.-', label="elem-max-norm")
    ax.loglog(Ms, err_tup[3], '.-', label="elem-1-norm")
    cmplxn(ax, 1, err_tup[2][-1])
    cmplxn(ax, 2, err_tup[0][-1])
    cmplxn(ax,3, err_tup[3][-1])
    ax.legend()
    ax.set_xlabel("rows or columns")
    fig.savefig(savedir / 'rounderr_matmul.png', format='png')

def plot_time_matmul(tM):
    fig, ax = plt.subplots()
    ax.loglog(Ms, tM, label="time [s]")    
    cmplxn(ax, 2, tM[-1])
    cmplxn(ax, 3, tM[-1])
    ax.legend()
    ax.set_xlabel("rows or columns")
    fig.savefig(savedir / 'cmplx_matmul.png', format='png')

def cmplxn(ax, n, scale=1):
    cmplx = [m**n/max(Ms)**n*scale*2 for m in Ms]
    ax.loglog(Ms, cmplx, '--', label="O(n^"+str(n)+")")

## Tests

# checking that it is not elementwise multiplication
def test_mat_multiply():
    I = np.asmatrix([[1, 0], [0, 1]])
    A = make_mat(2)
    IA = mat_multiply(I,A)
    assert IA[0, 1] == A[0, 1]
    assert IA[1, 0] == A[1, 0]

# cross-checking with vector instead of matrix
# found that error depends on the norm!
# found that rounding error of elementwise operation is O(n^0)
def test_vec_multiply():
    make_vec = lambda M: np.random.rand(M)
    # elementwise multiplication
    vec_multiply = lambda vec1, vec2: vec1 * vec2
    vec_round = lambda vec: np.asarray(vec, dtype=np.single)
    S = 10
    norm_arr = [
        lambda vec: np.sqrt(sum([v**2 for v in vec])),
        lambda vec: sum(abs(vec)),
        lambda vec: max(abs(vec)),
        lambda vec: sum([v**2 for v in vec]),
    ]
    err_tup = np.transpose(find_rounderr_sweep(vec_multiply, make_vec, vec_round, norm_arr, S, Ms))
    fig, ax = plt.subplots()
    ax.loglog(Ms, err_tup[0], '.-', label="2-norm")
    ax.loglog(Ms, err_tup[1], '.-', label="1-norm")
    ax.loglog(Ms, err_tup[2], 'o-', label="max-norm")
    ax.loglog(Ms, err_tup[2], '--', label="2-norm-square")
    cmplxn(ax, 0.25, err_tup[2][-1])
    cmplxn(ax, 0.5, err_tup[0][-1])
    cmplxn(ax, 1, err_tup[1][-1])
    ax.legend()
    ax.set_xlabel("rows or columns")
    fig.savefig(savedir / 'rounderr_vec.png', format='png')

if __name__ == "__main__":
    print("doing tests")
    test_mat_multiply()
    test_vec_multiply()
    print("finding rounding error")
    plot_matmul_rounderr(find_matmul_rounderr())
    print("finding complexity")
    plot_time_matmul(time_matmul())
    print("done! results saved to "+savedir.as_posix())