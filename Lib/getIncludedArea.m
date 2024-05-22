% Compute area between two curves (trajectories):
% A =  sum_N || x_demo_n - x(s_n) ||
% where the arc length s_n = ||x_demo_n - x_demo_(n-1)||
% N is the number of training samples
% 
% Author: Matteo Saveriano
% Date: 19.11.15

function [Area, realPtsInd] = getIncludedArea(xReal, xDemo)
    N = length(xDemo);
    T = length(xReal);
    
    % Start points are the same
    Area = norm(xDemo(:,1)-xReal(:,1));
    realPtsInd = zeros(1,N);
    realPtsInd(1) = 1;
    lastInd = 1;
    for trainIt=2:N        
        % Find next point in xReal
        arcLen = norm(xDemo(:,trainIt) - xDemo(:,trainIt-1));
        i = 1;
        while((lastInd+i)<T-2 && norm(xReal(:,lastInd+i) - xReal(:,lastInd))<arcLen)
            i = i+1;
        end
                
        if(lastInd+i>T)
            ee=1;
        end
        
        lastInd = min(lastInd + i,T);

        
        % Compute distance
        Area = Area + norm(xDemo(:,trainIt)-xReal(:,lastInd));
        
        % Store real points
        if(nargout > 1)
            realPtsInd(trainIt) = lastInd;
        end
    end
end