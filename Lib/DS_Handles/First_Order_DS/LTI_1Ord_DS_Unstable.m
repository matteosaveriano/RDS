function [xd, Vx] = LTI_1Ord_DS_Unstable(x, K)   
    % Unstable assuming K positive definite
    xd  = K*x;
    if(nargout>1)
        Vx = x; 
    end
end
