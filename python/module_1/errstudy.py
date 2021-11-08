import timeit

# Finds average round-off error over data sizes
# function test_operator(data1, data2)
# function test_data_maker(data_size)
# function test_operator_rounder(data)
# array of functions norm(data)
# int sample_size
# iterable data_sizes
# Returns [[size1 x norms...], [size2 x norms...], ...]
def find_rounderr_sweep(op, maker, op_round, norm_arr, S, Ms):
    err = [None for _ in range(len(Ms))]
    for i, M in enumerate(Ms):
        err[i] = find_rounderr(op, maker, op_round, norm_arr, S, M)
    return err

def find_rounderr(op, maker, op_round, norm_arr, S, M):
    Merr = [0 for _ in range(len(norm_arr))]
    for _ in range(S):
        data1 = maker(M); data2 = maker(M)
        outp = op(data1, data2)

        data1_lp = op_round(data1); data2_lp = op_round(data2)
        outp_lp = op(data1_lp, data2_lp)
        
        for j, norm in enumerate(norm_arr):
            Merr[j] += norm(outp_lp - outp) / S
    return Merr

# Times over data sizes
def time_sweep(op, maker, S, Ms):
    tM = [0 for _ in range(len(Ms))]
    for i, M in enumerate(Ms):
        data1 = maker(M); data2 = maker(M)
        # doc says: returns the time it takes to execute the main statement a number of times
        tM[i] = timeit.timeit(lambda: op(data1, data2), number=S) / S        
    return tM