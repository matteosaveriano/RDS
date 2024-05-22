% Multi Dimension Dynamic Time Warping Algorithm
% Dist is unnormalized distance between t and r
% D is the accumulated distance matrix
% k is the normalizing factor
% w is the optimal path
% t is the vector you are testing against
% r is the vector you are testing
% The orinal code by T. Felty [10]
% Modify by Parinya Sanguansat
function [Dist, D, k, w] = MD_DTW(t, r, getPath)

if(nargin < 3)
    getPath = 1;
end

[~,N]=size(t);
[rows,M]=size(r);

% This is the only difference between DTW and MD_DTW
d = 0;
for i=1:rows
    tt = t(i,:);
    rr = r(i,:);
    % Zero mean and unit variance for each dimension
%     std_t = std(tt);
%     if(std_t>1e-10)
%         tt = (tt-mean(tt))/std_t;
%     else
%         tt = (tt-mean(tt));
%     end
%     std_r = std(rr);
%     if(std_r>1e-10)
%         rr = (rr-mean(rr))/std_r;
%     else
%         rr = (rr-mean(rr));
%     end
    d = d + (repmat(tt(:),1,M) - repmat(rr(:)',N,1)).^2;
end

% The DTW algorithm
[Dist, D, k, w] = executeDtw(N, M, d, getPath);