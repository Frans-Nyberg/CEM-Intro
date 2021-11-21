function trips = make_trips(M,N,P)
% Makes (m,n,p) row vectors in a matrix

trips = [];
for m=0:M-1
    for n=0:N-1
        for p=0:P-1
            if ~(m==0 && n==0 && p==0)
                tr = [m,n,p];
                trips = [trips; tr];
            end
        end
    end
end

end