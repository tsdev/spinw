function varargout = sw_plotsf(sFact, varargin)
% fHandle = SW_PLOTSF(sFact, 'Option1', Value1, ...) plots the structure
% factor in the selected range.
%
% Options:
%
% range     Plotting range in r.l.u., dimensions are [3 2]. If the lower
%           and upper range is equal along one of the directions, the
%           signal will be integrated along that dimension. Default is 
%           [0 0;0 1;0 1].
% type      Type of plot.
%               'sf'    plot calculated structure factor. (default)
%               'perp'  plot the perpendicular component of the structure
%                       factor to Q.
%
% See also SW, SW.STRUCTFACT.
%

if nargin == 0
    help sw_plotsf
    return
end

inpForm.fname  = {'range'       'type' };
inpForm.defval = {[0 0;0 1;0 1] 'sf'   };
inpForm.size   = {[3 2]         [1 -1] };

param = sw_readparam(inpForm, varargin{:});
range = param.range;

% Dimensions to integrate along.
dInt   = find((range(:,1)-range(:,2))==0);
dRange = find((range(:,1)-range(:,2))~=0);

switch param.type
    case 'sf'
        F2 = sFact.int;
    case 'perp'
        F2 = sFact.perp;
    otherwise
        warning('sw:sw_plotsf:WrongParam','param.type is wrong, using ''sf''!');
        F2 = sFact.int;
end

for ii = 1:length(dInt)
    F2 = sum(F2,dInt(ii));
end

% TODO: Convolute with Gaussian function
if length(dRange) == 2
    F2 = squeeze(F2);
    fHandle = figure;
    pAxis   = {sFact.h sFact.k sFact.l};
    lAxis   = {'k_x (A^{-1})' 'k_y (A^{-1})' 'k_z (A^{-1})'};
    xAxis   = pAxis{dRange(1)}/2/pi;
    yAxis   = pAxis{dRange(2)}/2/pi;
    
    % the axis in reciprocal space (Angstrom^-1)
    Q1 = zeros(1,3);
    Q1(dRange(1)) = 1;
    Q1 = Q1/sFact.basisvector*2*pi;
    
    Q2 = zeros(1,3);
    Q2(dRange(2)) = 1;
    Q2 = Q2/sFact.basisvector*2*pi;
    
    Q3 = cross(Q1,Q2);
    
    % Mirror F2 matrix for all four quarters.
    F2a = rot90(F2(2:end,2:end),2);
    F2b = flipud(F2(2:end,:));
    F2c = fliplr(F2(:,2:end));
    F2d = F2;
    F2ext = [F2a F2b;F2c F2d];
    xAxis = [fliplr(-xAxis) xAxis(2:end)];
    yAxis = [fliplr(-yAxis) yAxis(2:end)];
    
    
    % rectangular grid in r.l.u. space
    [xxRlu, yyRlu] = ndgrid(xAxis,yAxis);
    % convert into real space (Angstrom^-1)
    xxA = xxRlu*Q1(1)+yyRlu*Q2(1);
    yyA = xxRlu*Q1(2)+yyRlu*Q2(2);
    zzA = xxRlu*Q1(3)+yyRlu*Q2(3);
    
    sf = surf(xxA,yyA,zzA,F2ext);
    set(sf,'EdgeAlpha',0);
    title('Magnetic structure factor')
    xlabel(lAxis{1});
    ylabel(lAxis{2});
    zlabel(lAxis{3});
    
    range(dInt,1) = -1;
    range(dInt,2) =  1;
    axis(reshape(range',[1,6]));
    Q3(Q3<1e-10) = 0;
    view(Q3);
    colorbar;
    cc = caxis;
    caxis([0 cc(2)]);
else
    error('sw:sw_plotsf:Dim','Only 2D plots are supported!');
end

if nargout>0
    varargout{1} = fHandle;
end
end