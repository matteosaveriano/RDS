function [xd, Vx] = LTI_2Ord_DS_goal(x, K, D, goal)
    [spaceDim, demoNum] = size(x);
    spaceDim = spaceDim/2;
       
    v = x(spaceDim+1:2*spaceDim,:);
    
    a   = zeros(spaceDim, demoNum);
    Vx  = zeros(2*spaceDim, demoNum); % Gradient of Lyaounov func. wrt x
    for i=1:demoNum
        a(:,i)  = K*(goal(1:spaceDim)-x(1:spaceDim,i)) - D*v(:,i);
        Vx(:,i) = [K*(x(1:spaceDim,i)-goal(1:spaceDim)); v(:,i)];
    end
    
    xd = [v; a];
end