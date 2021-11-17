
def iter_trips(M,N,P):
    for m in range(M):
        for n in range(N):
            for p in range(P):
                if m or n or p:
                    yield (m,n,p)
