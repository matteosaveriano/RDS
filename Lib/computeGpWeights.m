% This function computes the truncated GP weights (trunc((Kxx+sI)^(-1) * Kxs))
%
% Author: Matteo Saveriano
% Date: 02.11.15
%
% INPUT: GPR         -> Gaussian Process Regression structure
%        xStar       -> Current point
%        truncValues -> Parameters of the truncation function (see Eq. (11)) 
%
% OUTPUT: alpha -> Vector of truncated GP weights (trunc((Kxx+sI)^(-1) * Kxs))

function alpha = computeGpWeights(GPR, xStar, truncValues)
    % Compute training set covariance matrix
    Kxx = squaredExponentialCov(GPR);            
    
    % Compute test covariance
    Kxs = squaredExponentialCov(GPR, xStar);   
    
    % Compute weights (Kxx)^(-1) * Kxs
    alpha = Kxx\Kxs;
    
    % Truncate weights
    N = length(alpha);
    if(nargin < 3)
        a = 0.01;
        r = 0.01;
    else
        a = truncValues.alphaMin;
        r = truncValues.rho;
    end
    
    for i=1:N
        if(alpha(i) <= a)
            alpha(i) = 0;
        elseif(alpha(i) <= a+r)
            tr = 2*pi*(alpha(i)-a)/2*r;
            truncVal = 0.5*(1+sin(tr-pi/2));
            alpha(i) = truncVal * alpha(i);
        end
    end
end