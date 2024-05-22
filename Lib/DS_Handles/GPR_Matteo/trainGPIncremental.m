% This function train a 2D GP model incrementally
%
% Author: Matteo Saveriano
% Date: 02.11.15
%
% INPUT: GPR_old      -> Vector (2D) of initialized GP structures
%        trainingIn   -> Matrix (DxM) of the demonstrated state trajectory
%        trainingData -> A structure that contains:
%                         - scale -> Array (1xM) of scaling factors  
%                         - angle -> Array (1xM) of rotation angles
%        thresh       -> If(costFun(pt)>thresh) add pt to GP 
%
% OUTPUT: GPR -> Learned GP model (struct)

function GPR = trainGPIncremental(GPR_old, trainingIn, trainingData, threshVel, verbouse)
    if(nargin < 5)
        verbouse = false;
    end
    if(nargin < 4)
        threshVel = 0.1;
    end
    
    dataLen   = size(trainingIn, 2);
    numAddPts = 0;
    startIndx = 1;
    D         = length(GPR_old);
    
    % Save old points
    GPR = GPR_old;
    % Add at least one point in case
    if(isempty(GPR(1).data.in))
        for p=1:D
            GPR(p).data.in(:,1)  = trainingIn(:,1);
            GPR(p).data.out(1,1) = trainingData.deltaXd(p,1)';
        end
        % Pre-comute GP regression gains assuming all the GPR with same hyperparameters
        Kxx = squaredExponentialCov(GPR(1));  % compute training set covariance matrix

        for iter=1:D
            GPR(iter).alphaGain  = Kxx\GPR(iter).data.out;
        end
            
        numAddPts = 1;
        startIndx = 2;
    end
        
        % Iteratively add new points
    for i=startIndx:dataLen
        costPos = getCostPosition(GPR, trainingIn(:,i), trainingData.deltaXd(:,i));
        if((costPos > threshVel)) % Add new point
            for p=1:D
                GPR(p).data.in(:,end+1)  = trainingIn(:,i);
                GPR(p).data.out(end+1,1) = trainingData.deltaXd(p,i)';
            end
            
            % Pre-comute GP regression gains assuming all the GPR with same hyperparameters
            Kxx = squaredExponentialCov(GPR(1));  % compute training set covariance matrix

            for iter=1:D
                GPR(iter).alphaGain  = Kxx\GPR(iter).data.out;
            end

            numAddPts = numAddPts + 1;
        end
    end
    
    if(verbouse)
        disp(['Training finished. Added ' num2str(numAddPts) ' new points over ' num2str(dataLen)])
    end  
   
end


% Compute the value of the cost function
function costPos = getCostPosition(GPR, newIn, newOut)
    predictOut = MD_GP_Regression(GPR, newIn);
        
    costPos = norm(newOut - predictOut);
end