function [Dist, D, k, w] = executeDtw(N, M, d, getPath)
% Dynamic Time Warping Algorithm execution
% N -> length of the vector you are testing against (t)
% M -> length of the vector you are testing (r)
% d -> squared distance matrix between t and r. Changing d you can use the
% standard DTW or the MD-DTW. See: P. SANGUANSAT - "Multiple 
% Multidimensional Sequence Alignment Using Generalized Dynamic Time
% Warping" - WSEAS Transactions on Mathematics

D = zeros(size(d));
D(1,1) = d(1,1);

for n=2:N
    D(n,1) = d(n,1) + D(n-1,1);
end
for m=2:M
    D(1,m) = d(1,m) + D(1,m-1);
end

for n=2:N
    for m=2:M
        % this double 'min' construction improves in 10-fold the Speed-up.
        D(n,m) = d(n,m) + min(D(n-1,m), min(D(n-1,m-1), D(n,m-1)));
    end
end

Dist = D(N,M);
n = N;
m = M;
k = 1;
w = [];
if(getPath == 1)
    w(1,:) = [N,M];
    while((n+m) ~= 2)
        if((n-1) == 0)
            m = m-1;
        elseif((m-1) == 0)
            n = n-1;
        else
            [~,number] = min([D(n-1,m), D(n,m-1), D(n-1,m-1)]);
            switch number
                case 1
                    n = n-1;
                case 2
                    m = m-1;
                case 3
                    n = n-1;
                    m = m-1;
            end
        end
        k = k+1;
        w = cat(1,w,[n,m]);
    end
end