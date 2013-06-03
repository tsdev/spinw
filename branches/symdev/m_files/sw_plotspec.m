function [fHandle0, pHandle0] = sw_plotspec(spectra, varargin)
% [fHandle, pHandle] = SW_PLOTSPEC(spectra, 'option1', value1 ...) plots
% spin wave spectra.
%
% Options:
%
% mode      Choose the type of plot:
%               1   dispersion,
%               2   absolute value of the correlation functions,
%               3   convoluted spectra.
%               4   FANCY PLOT MODE, only axLim parameter can be defined.
%                   (default)
% mapplot   Whether to plot convoluted spectra using spec1d mapplot.
%           Default is false.
% imag      Whether to plot the imaginary values of the dispersion
%           and the correlation functions. For convoluted spectra, if true,
%           the imaginary part is plotted. Default is false.
% aHandle   Handle of the axis object for plotting, if undefined the
%           previous plot sw_plotspec plot window is used.
% colorbar  Plot colorbar for dispersion and intensity, default is true.
% nCol      Number of colors in the colormap, default is 500.
% dashed    Whether to plot dashed vertical line between multiple linear
%           scans. Defult is false.
% convE     FWHM value of convoluted Gaussian in energy to simulate finite
%           energy resolution. Only works for mode=3. If zero, no
%           convolution performed. Default is 0.
% fontSize  Font size on the plot, default is 14 pt.
% colormap  Colormap for plotting, default is @fireprint (good for color
%           and B&W printing) for single plot and for multiple plot it will
%           be a continuous scale from white to different color. This is
%           the 'auto' mode. Also colormap can be given directly using
%           standard colormaps, like @jet. To overplot multiple spectras
%           'colormap' option will be a matrix, with dimensions [3 nConv],
%           where every column defines a color for the maximum intensity.
%           It is also used for plotting dispersion curves. In case a
%           single color all dispersion curves have the same color (e.g.
%           [255 0 0] for red), or as many colors as dispersion curves
%           (dimensions are [3 nMode]), or any colormap can be given, like
%           @jet. In this case evry mode will have different colors, the
%           color is determined from the index of the mode. Default is
%           'auto'.
% sortMode  Sorting the modes before plotting. Default is false.
% axLim     Upper limit for y axis (mode 1,2) or z axis (mode 3), default
%           is 'auto'. For color plot of multiple cross section the c axis
%           cannot be changed after the plot.
% legend    Whether to plot legend for multiple convoluted spectras,
%           default is true.
% title     Whether to plot figure title, default is true.
% twin      Select which twins to plot for omega plots, default plots all
%           twins, dimensions are [1 nTwinToPlot].
% lineStyle Line style for line plots (dispersion and intensity), default
%           is {'-' 'o-' '--'}. For example '--' gives dashed lines.
% lineWidth Line width of line plots, default is 0.5 point.
%
% Output:
%
% fHandle   Handle of the plot figure.
% pHandle   Handle of the graphics objects on the figure.
%
% See also SW.PLOT, SW.SPINWAVE, SW.SWINC.
%

if nargin==0
    help sw_plotspec;
    return;
end

inpForm.fname  = {'mode' 'mapplot' 'imag' 'aHandle' 'colorbar' 'dashed' };
inpForm.defval = {4      false     false   0         true       false    };
inpForm.size   = {[1 1]  [1 1]     [1 1]  [1 1]     [1 1]      [1 1]    };

inpForm.fname  = [inpForm.fname  {'convE' 'fontSize' 'colormap' 'axLim'}];
inpForm.defval = [inpForm.defval {0       14         'auto'     'auto' }];
inpForm.size   = [inpForm.size   {[1 1]   [1 1]      [-1 -2]    [1 -3] }];

inpForm.fname  = [inpForm.fname  {'legend' 'title' 'nCol' 'twin'     }];
inpForm.defval = [inpForm.defval {true     true    500    zeros(1,0) }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]   [1 1]  [1 -4]     }];

inpForm.fname  = [inpForm.fname  {'lineStyle'     'lineWidth' 'sortMode'}];
inpForm.defval = [inpForm.defval {{'-' 'o-' '--'} 0.5         false     }];
inpForm.size   = [inpForm.size   {[1 -5]          [1 1]       [1 1]     }];
param = sw_readparam(inpForm, varargin{:});

% select twins for omega plot
param.twin = round(param.twin);
if isfield(spectra,'omega') && iscell(spectra.omega)
    nTwin      = numel(spectra.omega);
    if isempty(param.twin)
        param.twin = 1:nTwin;
    end
    param.twin = param.twin((param.twin<=nTwin) & (param.twin>0));
    if isempty(param.twin)
        warning('sw:sw_plotspec:WrongInput','Number of twins is wrong, plotting all twins!');
        param.twin = 1:nTwin;
    end
    nTwin      = numel(param.twin);
else
    nTwin = 1;
    param.twin = 1;
end

% select twins for convoluted plots
if iscell(spectra.swConv)
    % number of convoluted spectras to plot
    nTwinS      = size(spectra.swConv,2);
    param.twinS = param.twin((param.twin<=nTwinS) & (param.twin>0));
    if isempty(param.twinS) && (nTwinS>1)
        warning('sw:sw_plotspec:WrongInput','Number of twins is wrong, plotting all twins!');
        param.twinS = 1:nTwinS;
    end
    nTwinS      = numel(param.twinS);
else
    nTwinS      = 1;
    param.twinS = 1;
end

% Determine powder mode
powmode = false;
if numel(spectra.hklA)==length(spectra.hklA)
    powmode = true;
end

if param.mode == 4
    % PLOT EASY PEASY
    
    Eres = (spectra.Evect(end) - spectra.Evect(1))/50;
    [fHandle, pHandle1] = sw_plotspec(spectra,'ahandle',gca,'mode',3,'conve',Eres,...
        'dashed',true,'colorbar',false,'axLim',param.axLim);
    if ~powmode
        hold on
        [~, pHandle2] = sw_plotspec(spectra,'mode',1,'ahandle',gca,'colorbar',false,...
            'dashed',false,'title',false,'legend',false,'imag',false);
    end
    
    if nargout >0
        fHandle0 = fHandle;
    end
    if nargout>1
        pHandle0 = [pHandle1 pHandle2];
    end
    return
end

% Label of the x-axis
if powmode
    % powder mode
    xLabel  = 'Momentum transfer (A^-1)';
    xAxis   = spectra.hklA;
else
    [xLabel, xAxis] = sw_label(spectra);
end

yAxis  = spectra.Evect;
yLabel = 'Energy transfer (meV)';

if any(param.aHandle)
    axes(param.aHandle);
    fHandle = get(gca,'Parent');
else
    fHandle = sw_getfighandle('sw_spectra');
    if isempty(fHandle)
        fHandle = figure;
        set(fHandle,'Tag','sw_spectra');
    end
end

% Plotting styles for commensurate/incommensurate structures.
if iscell(param.lineStyle)
    if numel(param.lineStyle) ~= 3
        param.lineStyle = repmat(param.lineStyle(1),[1 3]);
    end
else
    param.lineStyle = repmat({param.lineStyle},[1 3]);
end


if ~powmode
    % Convert the convoluted intensities into cell array.
    if ~iscell(spectra.convmode)
        swInt    = {spectra.swInt};
        swConv   = {spectra.swConv};
        convmode = {spectra.convmode};
    else
        swInt    = spectra.swInt(:,param.twinS);
        swConv   = spectra.swConv(:,param.twinS);
        convmode = spectra.convmode;
    end
    % number of different convoluted cross sections
    nConv = numel(convmode);
    % number of convoluted plots
    nPlot  = nTwinS * nConv;
    
    % package all fields into cells for easy looping over twins
    if ~iscell(spectra.omega)
        omega = {spectra.omega};
    else
        omega = spectra.omega(:,param.twin);
    end
    
    nMode   = size(omega{1},1);
    nMagExt = spectra.obj.nmagext;
    % Defines colors for plotting modes.
    %colors  = flipud(fireprint(nMode+2));
    if isa(param.colormap,'function_handle')
        colors = flipud(param.colormap(nMode+2));
    else
        if strcmpi(param.colormap,'auto')
            colors = flipud(fireprint(nMode+2));
        else
            if numel(param.colormap) == 3
                param.colormap = param.colormap(:);
                colors = repmat(param.colormap',nMode+2,1)/255;
            elseif (size(param.colormap,1) == nMode) && (size(param.colormap,2)==3)
                colors = [0 0 0; param.colormap; 0 0 0]/255;
            elseif (size(param.colormap,2) == nMode) && (size(param.colormap,1)==3)
                colors = [0 0 0; param.colormap'; 0 0 0]/255;
            else
                error('sw:sw_plotspec:ColormapError','The dimensions of the colormap should be [3 nMode]');
            end
        end
    end
    colors  = colors(2:(end-1),:);
    
    modeList = nMode/(2*nMagExt);
    if modeList == 1
        if param.imag
            lLabel = {'Real' 'Imaginary'};
        else
            lLabel = {'Real'};
        end
    elseif modeList == 3
        if param.imag
            lLabel = {'Q+k_m' 'Q' 'Q-k_m' 'Imaginary'};
        else
            lLabel = {'Q+k_m' 'Q' 'Q-k_m'};
        end
    else
        error('sw:sw_plotspec:NumberOfModes','Wrong number of spin wave modes!');
    end
else
    nPlot    = 1;
    swConv   = {spectra.swConv};
    nConv    = 1;
    convmode = {spectra.convmode};
    param.legend = false;
end

if powmode && (param.mode~=3)
    warning('sw:sw_plotspec:PowMode','Powder spectra, only convoluted spectra can be plotted!');
    param.mode = 3;
end

hPlot = [];
hold on

switch param.mode
    case 1
        % Line plot of dispersion
        axis0 = [xAxis(1) xAxis(end) yAxis(1) yAxis(end)];
        titleStr0 = 'Spin wave dispersion: \omega(Q)';
        % loop over the twins
        for tt = 1:nTwin
            plotr = real(omega{1,tt});
            ploti = imag(omega{1,tt});
            if param.sortMode
                plotr = sort(plotr,1);
                ploti = sort(ploti,1);
            end
            % loop over all spin wave modes
            for ii = 1:nMode
                incIdx = ceil(ii/2/nMagExt);
                hPlot(end+1)    = plot3(xAxis,abs(plotr(ii,:)),xAxis*0+1e5,param.lineStyle{incIdx},...
                    'Color', colors(ii,:),'LineWidth',param.lineWidth); %#ok<*AGROW>
                hLegend(incIdx) = hPlot(end);
                if param.imag
                    hPlot(end+1)        = plot(xAxis,abs(ploti(ii,:)),'ro-');
                    hLegend(modeList+1) = hPlot(end);
                end
            end
        end
        
    case 2
        % Line plot of cross sections but only the first cell array element
        axis0 = [xAxis(1) xAxis(end) 0 1];
        yLabel = 'Intensity (arb. u.)';
        if param.imag
            titleStr0 = 'Intensity of the spin-spin correlation function: Im ';
        else
            titleStr0 = 'Intensity of the spin-spin correlation function: Re ';
        end
        % loop over the twins
        for tt = 1:nTwinS
            for jj = 1:nConv
                if param.imag
                    plotr = imag(swInt{jj,tt});
                else
                    plotr = real(swInt{jj,tt});
                end
                for ii = 1:nMode
                    hPlot(end+1) = plot3(xAxis,abs(plotr(ii,:)),xAxis*0+1e5,param.lineStyle{mod(jj-1,3)+1},...
                        'Color', colors(ii,:),'LineWidth',param.lineWidth); %#ok<*AGROW>
                    if ii == nMode
                        hLegend(jj) = hPlot(end);
                    end
                end
            end
            if param.imag
                
                for jj = 1:nConv
                    ploti = imag(swInt{jj,tt});
                    for ii = 1:nMode
                        hPlot(end+1) = plot3(xAxis,abs(ploti(ii,:)),xAxis*0+2e5,'ro-');
                        hLegend(nConv+1) = hPlot(end);
                    end
                end
                
            end
        end
end

if param.mode < 3
    if strcmpi(param.axLim,'auto')
        autAxis = axis;
        axis([sort(axis0(1:2)) autAxis(3:4)]);
    else
        axis([sort(axis0(1:2)) param.axLim]);
    end
    box on
end

if param.mode == 1
    if param.legend
        legend(hLegend,lLabel{:});
    end
    
    if param.colorbar
        cHandle = colorbar;
        set(get(cHandle,'ylabel'),'String', 'Index of modes');
        colormap(colors);
        caxis([1 nMode]);
    end
end

if param.mode == 3
    
    % filter out imaginary, inf and NaN values
    for ii = 1:nPlot
        if param.imag
            swConv{ii} = imag(swConv{ii});
        else
            swConv{ii} = real(swConv{ii});
        end
        swConv{ii}(isnan(swConv{ii})) = 0;
        swConv{ii}(isinf(swConv{ii})) = 0;
    end
    
    % make colormap from white to given color
    % by defining a cell of functions, with dimension of [1 nPlot]
    if ~isa(param.colormap,'function_handle')
        if strcmpi(param.colormap,'auto')
            % for 'auto' mode an equally space hue values are created for
            % use with multiple sectra plot
            if nPlot>1
                param.colormap = hsv2rgb([(1:nPlot)'/nPlot ones(nPlot,2)])';
            else
                param.colormap = {@fireprint};
            end
        end
        if ~iscell(param.colormap)
            if (size(param.colormap,1) ~= 3) || (size(param.colormap,2)<nPlot)
                error('sw:sw_plotspec:ColormapError','The dimensions of the colormap should be [3 nPlot]');
            end
            tHandle = cell(1,nPlot);
            for ii = 1:nPlot
                tHandle{ii} = @(numSteps)makecolormap(param.colormap(:,ii)'/255,[1 1 1],numSteps);
            end
            param.colormap = tHandle;
        end
    else
        param.colormap = repmat({param.colormap},[1 nPlot]);
    end
    
    % Gaussian energy convolution kern
    if param.convE>0
        sG = param.convE/2.35482;
        x0 = spectra.Evect;
        dx = (x0(end)-x0(1))/(length(x0)-1);
        nG = ceil(3*sG/dx);
        xG = (-nG*dx):dx:(nG*dx);
        % Gaussian normalised intensity
        fG = exp(-(xG/sG).^2/2);
        fG = fG/sum(fG);
    else
        fG = 1;
    end
    
    % c axis limit
    if strcmpi(param.axLim,'auto')
        % determine maximum intensity automatically
        zi = reshape(cell2mat(swConv),1,[]);
        
        zi(isnan(zi)) = [];
        posLim = sort([zi(zi>0) 0],'descend');
        negLim = sort([zi(zi<0) 0],'ascend');
        axLim  = [negLim(ceil(end*8e-3)) posLim(ceil(end*8e-3))];
        axMagn = 10.^(floor(log10(abs(axLim))));
        axLim  = ceil(axLim./axMagn).*axMagn;
        axLim(isnan(axLim)) = 0;
        cMaxMax = max(abs(zi));
        if axLim(2)-axLim(1) == 0
            axLim = [0 1];
        end
        
        if axLim(1) < 0
            axLim = [-max(abs(axLim)) max(abs(axLim))];
        end
%         zi(zi<=0) = [];
%         maxsort  = sort(zi,'descend');
%         maxsortS = maxsort(ceil(end*8e-3));
%         magni    = 10^(floor(log10(maxsortS)));
%         cMax     = ceil(maxsortS/magni)*magni;
%         cMaxMax  = maxsort(1);
%         
%         if cMax <= 0
%             cMax = 1;
%         end
    else
        axLim   = param.axLim;
        cMaxMax = max(abs(axLim));
        
        if numel(axLim) == 1
            axLim = [0 axLim];
        end
    end
    
    % convolute spectra with Gaussian to simulate finite energy resolution
    for ii = 1:nPlot
        iTemp = swConv{ii}';
        iTemp = conv2(iTemp,fG,'same');
        swConv{ii} = iTemp';
    end
    
    if nPlot == 1
        % single spectra
        imageDisp = swConv{1}';
        
        % Use surf to hide the NaN numbers
        [X, Y] = meshgrid(xAxis,yAxis);
        if cMaxMax <1e-6
            hSurf = surf(X,Y,imageDisp'*0);
            axLim = [0 1];
        else
            hSurf = surf(X,Y,imageDisp');
        end
        view(2);
        set(hSurf,'EdgeAlpha',0);
        colormap(flipud(param.colormap{1}(param.nCol)));
    else
        
        % multiple spectra
        vMat = zeros([size(swConv{1}),nPlot]);
        for ii = 1:nPlot
            vMat(:,:,ii) = swConv{ii};
        end
        cMat = sw_multicolor(vMat, param.colormap, axLim, param.nCol,true);
        % plot image piece-by-pice for the different Q directions
        if iscell(xLabel)
            xCut  = xLabel{end};
            nCut  = numel(xCut);
            hPlot = zeros(1,nCut);
            for ii = 2:nCut
                selIdx = xAxis>=xCut(ii-1) & xAxis<=xCut(ii);
                hPlot(ii-1) = image(xAxis(selIdx),yAxis,cMat(:,selIdx,:));
            end
            
        else
            hPlot = image(xAxis,yAxis,cMat);
        end
        set(gca,'YDir','normal');
    end
    
    caxis(axLim);
        
    if param.colorbar && (nPlot == 1)
        cHandle = colorbar;
        set(get(cHandle,'ylabel'),'String', 'Intensity (arb. u.)');
    end
    for ii = 1:nPlot
        lColor      = param.colormap{ii}(2);
        hLegend(ii) = plot(-1e5,-1e5,'Color',lColor(1,:));
    end
    
    axis0 = [xAxis(1) xAxis(end) yAxis(1) yAxis(end)];
    axis([sort(axis0(1:2)) axis0(3:4)])
    
    grid off
    set(gca,'Layer','top')
    titleStr0 = 'Convoluted ';
    if powmode
        titleStr0 = [titleStr0 'powder '];
    end
    titleStr0 = [titleStr0 'spectra: '];
end

ylabel(yLabel);

if param.mode > 1
    titleStr = cell(1,nConv);
    for ii = 1:nConv
        titleStr{ii} = sw_titlestr(convmode{ii});
        if param.imag
            titleStr{ii} = ['Im ' titleStr{ii}];
        else
            titleStr{ii} = ['Re ' titleStr{ii}];
        end
    end
    
    for ii = 1:nConv-1
        titleStr0 = [titleStr0 titleStr{ii} ', '];
    end
    titleStr0 = [titleStr0 titleStr{end}];
    
    if param.legend
        % put the names of the twins into the plot
        titleStr = repmat(titleStr,[1 nTwinS]);
        if nTwinS > 1
            for ii = 1:nPlot
                titleStr{ii} = ['Twin #' num2str(ceil(ii/nConv)) ' ' titleStr{ii}];
            end
        end
        % put label for imaginary lineplot
        if param.imag && param.mode==2
            titleStr = [titleStr(1:nConv) lLabel(2)];
        end
        warn_state = warning;
        warning('off','MATLAB:legend:IgnoringExtraEntries')
        legend(hLegend,titleStr{:},'FontSize',param.fontSize);
        warning(warn_state);
    end
end

if param.title
    title(titleStr0,'FontSize',param.fontSize);
end
drawnow;
hold off

if iscell(xLabel)
    xTickLoc = xLabel{end};
    set(gca,'XTick',xTickLoc)
    set(gca,'XTickLabel',xLabel(1:end-1));
    if param.dashed
        autAxis = axis;
        for jj=2:length(xLabel)-1
            line([1 1]*xTickLoc(jj),autAxis(3:4),'LineStyle','--','color',[0 0 0]);
        end
    end
    xlabel('Momentum (r.l.u.)');
else
    xlabel(xLabel);
end

if nargout == 1
    fHandle0 = fHandle;
elseif nargout == 2
    fHandle0 = fHandle;
    pHandle0 = hPlot;
end
end

function titleStr = sw_titlestr(convmode)
% creates plot title string from param.convmode string

titleStr = [convmode '}(\omega,Q)'];
titleStr = strrep(titleStr,'perp','\perp ');
titleStr = strrep(titleStr,'S','S^{');
titleStr = strrep(titleStr,'P','P^{');
titleStr = strrep(titleStr,'M','M^{');
titleStr = strrep(titleStr,'+','}(\omega,Q) +');
titleStr = strrep(titleStr,'-','}(\omega,Q) -');

end