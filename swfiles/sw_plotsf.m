function varargout = sw_plotsf(sFact, varargin)
% plots the structure factor in the selected Q range in 1D or 2D
%
% fHandle = SW_PLOTSF(sFact, 'Option1', Value1, ...)
%
% Options:
%
% range     Data range in inverse Angstrom, dimensions are [1 2] or [2 2]
%           for 1D and 2D plots respectively, default is the full data
%           range.
% log       Plot 10 based logarithmic intensity, default is false.
% colorbar  To show colorbar, default is true.
% plotStyle Options to the plot command, default is '-'.
%
% See also SW, SW.STRUCTFACT, SW_INTSF.
%

if nargin == 0
    help sw_plotsf
    return
end

inpForm.fname  = {'range' 'log' 'colorbar' 'plotStyle'};
inpForm.defval = {[]      false true       '-'        };
inpForm.size   = {[-1 2]  [1 1] [1 1]      [1 -2]     };
inpForm.soft   = {1       0     0          0          };

param = sw_readparam(inpForm, varargin{:});


if numel(sFact.axis)==3
    % single number
elseif numel(sFact.axis) == 2 
    % 1D plot
    fHandle = figure;
    sFact.axis
    plot(squeeze(sFact.int),'o-');
    
    
elseif sFact.axis == 0
    % 2Theta plot
    fHandle = figure;
    if param.log
        semilogy(sFact.hklAint,sFact.int,param.plotStyle);
        ylabel('Logarithm of Intensity (arb. units.)');
        title('Powder average of magnetic structure factor (LOG)')
    else
        plot(sFact.hklAint,sFact.int,param.plotStyle);
        ylabel('Intensity (arb. units.)');
        title('Powder average of magnetic structure factor')
    end
    
    xlabel('Momentum Transfer (A^{-1})');    
    
elseif  numel(sFact.axis) == 1
    % 2D plot
    % TODO: Convolute with Gaussian function
    
    pAxis = find(1:3~=sFact.axis);
    iPlot   = squeeze(sFact.int);
    fHandle = figure;
    labels  = {'k_x (A^{-1})' 'k_y (A^{-1})' 'k_z (A^{-1})'};
    xAxis   = sFact.hkl{pAxis(1)};
    yAxis   = sFact.hkl{pAxis(2)};
    
    % the axis in reciprocal space (Angstrom^-1)
    Q1 = zeros(1,3);
    Q1(pAxis(1)) = 1;
    Q1 = Q1/sFact.obj.basisvector*2*pi;
    
    Q2 = zeros(1,3);
    Q2(pAxis(2)) = 1;
    Q2 = Q2/sFact.obj.basisvector*2*pi;
    
    Q3 = cross(Q1,Q2);
    
    % Mirror F2 matrix for all four quarters.
    F2a = rot90(iPlot(2:end,2:end),2);
    F2b = flipud(iPlot(2:end,:));
    F2c = fliplr(iPlot(:,2:end));
    F2d = iPlot;
    F2ext = [F2a F2b;F2c F2d];
    xAxis = [fliplr(-xAxis) xAxis(2:end)];
    yAxis = [fliplr(-yAxis) yAxis(2:end)];
    
    
    % rectangular grid in r.l.u. space
    [xxRlu, yyRlu] = ndgrid(xAxis,yAxis);
    % convert into real space (Angstrom^-1)
    xxA = xxRlu*Q1(1)+yyRlu*Q2(1);
    yyA = xxRlu*Q1(2)+yyRlu*Q2(2);
    zzA = xxRlu*Q1(3)+yyRlu*Q2(3);
    
    if param.log
        F2ext = log10(F2ext+1e-10);
    end
    
    sf = surf(xxA,yyA,zzA,F2ext);
    set(sf,'EdgeAlpha',0);
    
    if param.log
        titleStr = 'Logarithm of the magnetic structure factor';
    else
        titleStr = 'Magnetic structure factor';
    end
    
    title(titleStr);
    xlabel(labels{1});
    ylabel(labels{2});
    zlabel(labels{3});
    
    pRange{1} = [min(xxA(:)) max(xxA(:))];
    pRange{2} = [min(yyA(:)) max(yyA(:))];
    pRange{3} = [min(zzA(:)) max(zzA(:))];
    pRange{sFact.axis} = [-1 1];
    
    if ~isempty(param.range)
        pRange{pAxis(1)} = param.range(1,:);
        pRange{pAxis(2)} = param.range(2,:);
    end
    axis([pRange{:}]);
    
    Q3(abs(Q3)<1e-10) = 0;
    view(Q3);
    cc = caxis;
    if cc(1) > 0
        caxis([0 cc(2)]);
    else
        caxis(cc);
    end
    
    if param.colorbar
        cHandle = colorbar;
        cLabelU = '(arb. units)';
        if param.log
            cLabel = ['log Intensity ' cLabelU];
        else
            cLabel = ['Intensity ' cLabelU];
        end
        
        set(get(cHandle,'ylabel'),'String',cLabel);
    end
    
    
end


if nargout>0
    varargout{1} = fHandle;
end
end