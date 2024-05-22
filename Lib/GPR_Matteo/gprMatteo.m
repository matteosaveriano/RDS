%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     INCREMENTAL LEARNING DS                             %
%                    Copyright 2015 - Matteo Saveriano                    %
%                                                                         %
%      This program is distributed under the terms of the GNU License     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gaussian Process Regression

function [mu, Sigma] = gprMatteo(GPR, xStar)  
    Kxx = squaredExponentialCov(GPR);  % compute training set covariance matrix

%     L = chol(Kxx)';             % cholesky factorization of the covariance
%     alpha = solve_chol(L', GPR.data.out);

    alpha = Kxx\GPR.data.out;

    [Kxs, Kss] = squaredExponentialCov(GPR, xStar);     %  test covariances

    mu = Kxs' * alpha;                                      % predicted means

    if nargout == 2
        v = L\Kxs;
        Sigma = Kss - sum(v.*v)'; % predicted covariances
    end  
end
