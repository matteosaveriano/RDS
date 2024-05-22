function [Dist, D, k, w] = DTW(t, r, getPath)
% Dynamic Time Warping Algorithm
% Dist    -> is unnormalized distance between t and r
% D       -> the accumulated distance matrix
% k       -> the normalizing factor
% w       -> the optimal path
% t       -> the vector you are testing against
% r       -> the vector you are testing
% getPath -> 1 - compute w, 0 - skip this step
% The orinal code by T. Felty [10]

if(nargin < 3)
    getPath = 1;
end

[~,N] = size(t);
[~,M] = size(r);

% This is the only difference between DTW and MD_DTW
d = (repmat(t(:),1,M)-repmat(r(:)',N,1)).^2;
%this replaces the nested for loops from above Thanks Georg Schmitz

% The DTW algorithm
[Dist, D, k, w] = executeDtw(N, M, d, getPath);

%Dist = sqrt(Dist);

