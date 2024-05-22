function pos = reshapeDS(DSHandle, GPR, x0, simOptions, truncValues)
    for demoIt=1:size(x0,2)
        it = 1;
        xCurrMod = x0(:,demoIt);
        itM   = 0;
        while(it<simOptions.iterNum && norm(xCurrMod-simOptions.goal)>simOptions.tol)
            % Store data for plotting
            out_.xM(:,it) = xCurrMod;

            % Compute GP data
%             alpha = computeGpWeights(GPR(1), xCurrMod, truncValues);
%             w     = min(alpha' * GPR(1).data.out,1);           
            
            alpha     = computeGpWeights(GPR(2), xCurrMod, truncValues);
            par.xd(1,1) = alpha' * GPR(2).data.out;
            
            alpha     = computeGpWeights(GPR(3), xCurrMod, truncValues);
            par.xd(2,1) = alpha' * GPR(3).data.out;
            
            w = 1;
            if(it>simOptions.Tmax || norm(par.xd)<0.1 || norm(xCurrMod-simOptions.goal)<3)
                w = exp(-13*itM*simOptions.dt);
                itM = itM + 1;
            end
            
            
%             if(norm(par.xd)>0.01 && norm(xCurrMod-simOptions.goal)>3)
%                 w=1;
%             else
%                 w=0;
%             end

            % Compute next modulated position
            vCurrMod = DSHandle(xCurrMod);
            vCurrMod = vCurrMod + w*(par.xd - vCurrMod);
            xCurrMod = xCurrMod + vCurrMod * simOptions.dt;
%             xCurrMod = xCurrMod + w * par.xd * simOptions.dt;


            it = it + 1;
        end
        
        if(simOptions.plotResult)
            figure(1);
            hold on;
            plot(out_.xM(1,1:1:end), out_.xM(2,1:1:end), 'k','Linewidth', 3)
            plot(out_.xM(1,end), out_.xM(2,end), 'k.','Linewidth', 2,'MarkerSize', 50)
        end
        
        pos{demoIt} = out_.xM;
    end
end