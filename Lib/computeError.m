% Compute approximation error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Copyright (c) 2015 Matteo Saveriano, DHRI Lab, TUM,   %%%
%%%          80333 Munich, Germany, http://www.hri.ei.tum.de/en/home/   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The program is free for non-commercial academic use. Please contact the
% author if you are interested in using the software for commercial purposes.
% The software must not be modified or distributed without prior permission
% of the authors. Please acknowledge the authors in any academic publications
% that have made use of this code or part of it. Please use this BibTex
% reference:
% 

function error = computeError(DSHandle, GPR, xDemo, dt)
truncValues.alphaMin = 0.1;
truncValues.rho      = 0.1;

alpha       = computeGpWeights(GPR(1), xDemo', truncValues);
w(1,:)      = min(alpha' * GPR(1).data.out,1);
w(2,:)      = w(1,:); 
par.xd(1,:) = alpha' * GPR(2).data.out;
par.xd(2,:) = alpha' * GPR(3).data.out;

% Compute next modulated position
xd = DSHandle(xDemo); %compute outputs
[d, demoLen] = size(xd);
sumX0 = xDemo(:,1);
x = [];
for i=1:demoLen
    if(norm(xDemo(:,i))<2.)
        w(:,i) = [0;0];
    end
    xd(:,i) = xd(:,i) + w(:,i).*(par.xd(:,i) - xd(:,i));
    
    x = [x sumX0+xd(:,i)*dt];
end

for i=1:d 
    tmpError(i,:) = norm(x(i,:) - xDemo(i,:));
end

error = sum(tmpError);