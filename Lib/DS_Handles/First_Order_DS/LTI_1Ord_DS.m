function [xd, Vx] = LTI_1Ord_DS(x, K)   
    xd  = -K*x;
    if(nargout>1)
        Vx = x; % Assuming K positive definite
    end
end