%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     INCREMENTAL LEARNING DS                             %
%                    Copyright 2015 - Matteo Saveriano                    %
%                                                                         %
%      This program is distributed under the terms of the GNU License     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Multi-Dimensional Gaussian Process Regression

function mu = MD_GP_Regression(GPR, xStar)
    D  = length(GPR);
    N  = size(xStar, 2);
    mu = zeros(D,N);
    
    for j=1:N
        %  test covariances (assuming same training points for all GPR)
        Kxs = squaredExponentialCov(GPR(1), xStar(:,j)); 
        for i=1:D
            mu(i,j) = Kxs' * GPR(i).alphaGain;            % predicted means
        end
    end
end
