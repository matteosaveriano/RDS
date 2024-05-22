%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     INCREMENTAL LEARNING DS                             %
%                    Copyright 2015 - Matteo Saveriano                    %
%                                                                         %
%      This program is distributed under the terms of the GNU License     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the squared exponential covariance matrix

function [K1, K2] = squaredExponentialCov(GPR, xStar)
    lenScale = exp(GPR.logtheta(1));   % characteristic length scale
    sigmaF2  = exp(2*GPR.logtheta(2)); % signal squared variance
    sigmaN2  = exp(2*GPR.logtheta(3)); % noise squared variance
    
    nData = size(GPR.data.in', 2);
       
    if(nargin == 1) % compute training set covariance matrix (K1 = Kxx + sigmaN2*I)
        SDxx = sq_dist(GPR.data.in') / (2*lenScale);
        K1   = sigmaF2*exp(-SDxx) + sigmaN2*eye(nData);
    else % compute test set covariance
        % K1 = Kxs
        SDxs = sq_dist(GPR.data.in', xStar') / (2*lenScale);
        K1   = sigmaF2*exp(-SDxs);
        if(nargout == 2)
            K2 = sigmaF2*ones(size(xStar,1),1) + sigmaN2; % K2 = Kss + sigmaN2*I
        end
    end
end