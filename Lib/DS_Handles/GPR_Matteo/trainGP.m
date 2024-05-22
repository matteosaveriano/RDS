%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     INCREMENTAL LEARNING DS                             %
%                    Copyright 2015 - Matteo Saveriano                    %
%                                                                         %
%      This program is distributed under the terms of the GNU License     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Train GP

% inData -> a MxD matrix, where M is the number of samples, D the space
% dimension
%
% inData -> a MxD matrix, where M is the number of samples, D the space
% dimension

function GPR = trainGP(GPR, inData, outData)  
    D = length(GPR);
    % Add input and output data to the GPR struct
    for i=1:D
        GPR(i).data.in  = inData;
        GPR(i).data.out = outData(i,:);
    end
    
    % Pre-comute GP regression gains assuming all the GPR with same hyperparameters
    Kxx = squaredExponentialCov(GPR(1));  % compute training set covariance matrix
    
    for i=1:D
        GPR(i).alphaGain  = Kxx\GPR(i).data.out';
    end
end
