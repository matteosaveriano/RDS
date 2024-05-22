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

function GPR = trainGP(GPR_old, trainingIn, trainingData, threshPos, verbouse)
    if(nargin < 5)
        verbouse = false;
    end
    if(nargin < 4)
        threshPos = 0.1;
    end
    
    dataLen   = size(trainingIn, 2);
    numAddPts = 0;
    startIndx = 1;
    
    % Save old points
    GPR = GPR_old;
    % Add at least one point in case
    if(isempty(GPR(1).data.in))
        GPR(1).data.in(1,:) = trainingIn(:,1)';
        GPR(1).data.out(1,1)  = trainingData.act(1);
        numAddPts = 1;
        startIndx = 2;
    end
    if(isempty(GPR(2).data.in))
        for p=1:length(GPR)-1
            GPR(p+1).data.in(1,:)  = trainingIn(:,1)';
            GPR(p+1).data.out(1,1) = trainingData.xd(p,1)';
        end
        numAddPts = 1;
        startIndx = 2;
    end
        
        % Iteratively add new points
    for i=startIndx:dataLen
        costPos = getCostPosition(GPR(2:end), trainingIn(:,i)', trainingData.xd(:,i)');
        if((costPos > threshPos)) % Add new point
                GPR(1).data.in(end+1,:)  = trainingIn(:,i)';
                GPR(1).data.out(end+1,1) = trainingData.act(i);

                for p=1:length(GPR)-1
                    GPR(p+1).data.in(end+1,:)  = trainingIn(:,i)';
                    GPR(p+1).data.out(end+1,1) = trainingData.xd(p,i)';
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
    dim_ = length(GPR);
    predictOut = zeros(1, dim_);
    for i=1:length(GPR)
        predictOut(i) = gprMatteo(GPR(i), newIn);
    end
    
    costPos = norm(newOut - predictOut);
end