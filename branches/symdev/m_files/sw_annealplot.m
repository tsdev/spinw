function sw_annealplot(T, E, rate, param,fid)
% SW_ANNEALPLOT(T, E, rate, param,fid) displays information about the annealing
% procedure.
%

if nargin == 0
    help sw_annealplot;
    return
end

if  (param.verbosity > 0) && ~any(any(E))
    tic;
    dateStr = date;
    clockV  = clock;
    timeStr = sprintf('%02d:%02d:%02d',floor(clockV(4:6)));
    fprintf0(fid,'Starting time: %s %s.\n',dateStr,timeStr);
end

if (param.verbosity > 0) && ~isempty(E)
    MCTime = toc; tic;
    % Display current temperature, system energy, acceptance rate and number of
    % steps per temperature per spin.
    if size(E,2) == 1
        fprintf0(fid,'  T = %11.6f, E = %10.5f per moment, acc\\rej=%10.5f, Steps=%6d, Time=%10.3f s\n',T(end),...
            E(end),sum(rate)/length(rate),param.nMC,MCTime);
    else
        % Output for parallel tempering.
        fprintf0(fid,'  Time=%7.0f s\n',MCTime*param.nCycle);
    end
    
    if param.verbosity > 1
        % Handle of figure window to plot.
        hFigure = sw_getfighandle('sw_anneal');
        if isempty(hFigure)
            sw_annealfigure;
        end
        
        if size(E,2) == 1
            % Plot for simulated annealing.
            if ~isempty(E)
                nPlot = 100;
                ratePlot = zeros(1,nPlot);
                
                for ii = 1:nPlot
                    idx = ceil((ii-1)*length(rate)/nPlot+1):ceil(ii*length(rate)/nPlot);
                    ratePlot(ii) = sum(rate(idx))/length(idx)*100;
                end
                
                % Plots spin flip rate as a function of time.
                subplot(2,1,1); box on; grid on; hold on; pcolor = jet(1000);
                plot(linspace(0,100,nPlot),ratePlot,'Color',pcolor(randi(600,1)+200,:));
                title('Acceptance rate at the last temperature during thermalisation')
                xlabel('Time (percent)'); ylabel('Acceptance rate (percent)')
                axis([0 100 0 min(max(ratePlot)*3,100)+1e-8]);
                
                % Plots system energy as a function of temperature.
                subplot(2,1,2); hold off; semilogx(T, E, 'o-'); grid on;
                xlabel('Temperature (K)'); ylabel('Energy per extended magnetic unit cell (meV)');
            end
        elseif any(any(E))
            % Plot for parallel tempering.
            % Bin the energy distribution.
            E    = E(any(E'),:);
            nT   = length(T);
            minE = min(min(E));
            maxE = max(max(E));
            nX   = nT*15;
            x    = linspace(minE,maxE,nX);
            binE = zeros(nT,nX);
            
            for ii = 1:size(E,1)
                for jj = 1:nT
                    idx = floor((E(ii,jj)-minE)/(maxE-minE)*(nX-1))+1;
                    binE(jj,idx) = binE(jj,idx) + 1;
                end
            end
            binE = binE/max(max(binE));
            
            % Create random colors for the Gaussian curves.
            cList = jet(nT);
            cList = cList(mod((1:nT)*103,nT)+1,:);
            
            % Calculate overlaps.
            oLap = zeros(1,nT-1);
            for ii = 1:nT-1
                oLap(ii) = 2*sum(min([binE(ii,:); binE(ii+1,:)]))/sum(sum(binE(ii+(0:1),:)));
            end
            
            % Plot the energy distribution and overlaps.
            subplot(2,1,1)
            cla;
            for ii = 1:nT
                plot(x,binE(ii,:),'Color',cList(ii,:));
                hold on
            end
            pOL = plot(linspace(minE,maxE,nT-1),oLap,'ko-','LineWidth',1.5);
            legend(pOL,'Energy overlap of neighbouring temperatures')
            axis([minE maxE 0 1]);
            title(['Cycle #' num2str(size(E,1))]);
            xlabel('Energy per spin')
            ylabel('Normalized distribution')
            
            % Plots system energy as a function of temperature.
            subplot(2,1,2); hold off; semilogx(T, E(end,:), 'o-'); grid on;
            xlabel('Temperature (K)'); ylabel('Energy per extended magnetic unit cell (meV)');
            axis([min(T) max(T) minE maxE]);
        end
        drawnow;
    end
    if (size(E,2) == 1) && (T(end)<param.endT)
        dateStr = date;
        clockV  = clock;
        timeStr = sprintf('%02d:%02d:%02d',floor(clockV(4:6)));
        fprintf0(fid,'End time: %s %s.\n',dateStr,timeStr);
    end
    
end

end