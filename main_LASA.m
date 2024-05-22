% Reshaped Dynamical Systems
% 
% Transform a Linear DS into a complex motion
%
% Author: Matteo Saveriano
% Date: 22.05.24

close all
clear

addpath(genpath('./'));

%% Load demonstrations
motionNum = 30;
[demos, dt] = load_LASA_models('LASA_hand_writing_dataset/', motionNum);

% DS_handle = @(x) GMR(Priors,Mu,Sigma,x,1:d,d+1:2*d);
DS_handle = @(x) LTI_1Ord_DS(x, 3);

%% Parameters definition
% Simulation parameters
samplingRate          = 10;
simOptions.dt         = dt*samplingRate;
simOptions.iterNum    = 1000;
simOptions.tol        = 1e-3;
simOptions.dsOrder    = 1;

%% Incremental Learning
% Buil GPR structures
% GPR(1) window scale
GPR(1).logtheta = [log(2.); log(sqrt(1.0)); log(sqrt(.3))];
GPR(1).covfunc = {'covSum', {'covSEiso','covNoise'}};
GPR(1).data.in = [];
GPR(1).data.out = [];

% GPR(2)...GPR(N) desired state
GPR(2) = GPR(1);
GPR(3) = GPR(2);

demoNum = [1 3 5 7];
demoSize = size(demos{1}.pos, 2);
sInd = 1:samplingRate:demoSize;
figure(1)
hold on;
xm = []; xM = []; ym = []; yM = [];
for demoIt=demoNum   
    % Compute Learning Parameters
    xDemo  = demos{demoIt}.pos(:,sInd);
    xdDemo = demos{demoIt}.vel(:,sInd);
     
    % Plot demos
    plot(xDemo(1,:), xDemo(2,:), 'o','color',[1, 0, 1],'Linewidth', 2, 'MarkerSize', 5)
    
    trainingData.act = 1.05*ones(length(xDemo), 1);
    trainingData.xd  = xdDemo;

    % Learn GPs
    threshVel = 0.1;
    GPR = trainGP(GPR, xDemo, trainingData, threshVel, true);
    
    % Find figure box size
    xm = [xm min(demos{demoIt}.pos(1,sInd))];
    xM = [xM max(demos{demoIt}.pos(1,sInd))];
    ym = [ym min(demos{demoIt}.pos(2,sInd))];
    yM = [yM max(demos{demoIt}.pos(2,sInd))];
end

deltaBox = 5;
xmin = min(xm) - deltaBox;
xmax = max(xM) + deltaBox;
ymin = min(ym) - deltaBox;
ymax = max(yM) + deltaBox;

goal = [0; 0];
truncValues.alphaMin = 0.01;
truncValues.rho      = 0.01;
Tmax  = length(sInd) -5;

GPR_plot = GPR;
GPR_plot(1).logtheta = [log(30); log(sqrt(1.0)); log(sqrt(2.))];
GPR_plot(2).logtheta = [log(30); log(sqrt(1.0)); log(sqrt(2.))];
plotDSStreamLines(DS_handle, GPR_plot, [xmin,xmax,ymin,ymax], 'low')

%% Motion Reproduction
repIter = 1;
for demoIt=demoNum
    it = 1;
    x0 = demos{demoIt}.pos(:,1);
    xCurr = x0;
    vCurr = demos{demoIt}.vel(:,1);
    xCurrMod = x0;
    vCurrMod = vCurr;
    itM   = 0;
    while(it<simOptions.iterNum && norm(xCurrMod-goal)>simOptions.tol)
        % Store data for plotting
        out_.x(:,it) = xCurr;
        out_.xd(:,it) = vCurr;
        out_.xM(:,it) = xCurrMod;
        out_.xdM(:,it) = vCurrMod;

        % Compute next original position
        vCurr = DS_handle(xCurr);
        xCurr = xCurr + vCurr * simOptions.dt;

        % Compute GP data
        alpha = computeGpWeights(GPR(1), xCurrMod', truncValues);
        w     = min(alpha' * GPR(1).data.out,1);
        if(it>Tmax)
            w = exp(-13*itM*simOptions.dt);
            itM = itM + 1;
        end

%         alpha     = computeGpWeights(GPR(2), xCurrMod', truncValues);
        par.xd(1,1) = alpha' * GPR(2).data.out;
        par.xd(2,1) = alpha' * GPR(3).data.out;

        % Compute next modulated position
        vCurrMod = DS_handle(xCurrMod);
        vCurrMod = vCurrMod + w*(par.xd - vCurrMod);
        xCurrMod = xCurrMod + vCurrMod * simOptions.dt;

        % Store data for plotting
        out_.w(it) = w;
        out_.xGP(:,it) = par.xd;

        it = it + 1;
    end
    
    % Plot results
    figure(1);
    hold on;
%     plot(out_.x(1,:), out_.x(2,:), 'b','Linewidth', 2)
    plot(out_.xM(1,1:1:end), out_.xM(2,1:1:end), 'k','Linewidth', 3)
    %plot(out_.x(1,1), out_.x(2,1), 'k.','Linewidth', 2,'MarkerSize', 30)
       
    plot(out_.x(1,end), out_.x(2,end), 'k.','Linewidth', 2,'MarkerSize', 50)
    %title('Position [m]')
    
    % DTW distance
    for dim=1:2
        dist(dim) = DTW(out_.xM(dim,:), demos{demoIt}.pos(dim,:), 0);
    end
    totDist(repIter) = sum(dist);%/size(out_.xM,2); % /N to make points independent
    repIter = repIter + 1;
end

disp(['Total DTW Distance - Mean: ' num2str(mean(totDist)) ' - Std: ' num2str(std(totDist))])