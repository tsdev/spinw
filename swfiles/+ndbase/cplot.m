function hPlot = cplot(dat, varargin)
% produces interplolated color plot of scattered data
%
% hPlot = NDBASE.CPLOT(dat, 'option1', value1, ...)
%
% Input:
%
% dat       Data, either array of spec1d type or struct with fields x, y
%           and z.
%
% Options:
%
% x         X-values for spec1d data.
% scale     Scalar, determines the scaling between the x- and y-axis for
%           the interpolation. Default is 'auto'.
% npix      Number of pixels either scalar or 2 element vector
%           [npixx npixy].
% plot      Logical, if true a plot will be produced.
% scatter   If true, dots are plotted at the position of the given data
%           points. Default is true.
% convhull  If true, the color map is shown only between the given data
%           points (within the convex hull of the points), this is achieved
%           by giving nan value for all interpolated points outside of the
%           convex hull. Default is true.
% log       Plot natural logarithm of the values.
%
% Output
%
% hPlot     Handle of the plot object.
%

inpForm.fname  = {'x'    'scale' 'npix' 'plot' 'scatter' 'convhull' 'log'};
inpForm.defval = {[]     'auto'  200    true   true      true       false};
inpForm.size   = {[1 -2] [1 -3]  [1 -4] [1 1]  [1 1]     [1 1]      [1 1]};
inpForm.soft   = {true   false   false  false  false     false      false};

param = sw_readparam(inpForm, varargin{:});

if isa(dat,'spec1d')
    % read spec1d data + x values
    if numel(param.x)~=numel(dat)
        error('cplot:WrongInput','For spec1d data, option x needs to have nDat number of elements!');
    end
    
    x = [];
    y = [];
    z = [];
    
    for ii = 1:numel(dat)
        dd = struct(dat(ii));
        x  = [x;repmat(param.x(ii),[numel(dd.x) 1])];   %#ok<AGROW>
        y  = [y;dd.x];                                  %#ok<AGROW>
        z  = [z;dd.y];                                  %#ok<AGROW>
    end
    
else
    % take dat as structure
    x = dat.x(:);
    y = dat.y(:);
    z = dat.z(:);
end

if numel(param.npix) == 1
    npix = [1 1]*param.npix;
else
    npix = param.npix(1:2);
end

xi = linspace(min(x),max(x),npix(1));
yi = linspace(min(y),max(y),npix(2));

% scale
if ischar(param.scale) && strcmp(param.scale,'auto')
    scale = (max(y)-min(y))/(max(x)-min(x));
else
    scale = param.scale;
end

% do the interpolation
F = scatteredInterpolant(x(:)*scale,y(:),z(:),'natural');
[xi, yi] = meshgrid(xi,yi);
zi = F(xi*scale,yi);

% check if points are within the convex hull
if param.convhull
    hullIdx = convhull(x(:),y(:));
    zi(~inpolygon(xi(:),yi(:),x(hullIdx),y(hullIdx))) = nan;
end

if param.log
    zi = real(log(zi));
    zi(isnan(zi)) = 0;
end

% plot
if param.plot
    hPlot = pcolor(gca,xi,yi,zi);
    set(hPlot,'LineStyle','None')
    shading flat
    axis([min(x) max(x) min(y) max(y)]);
end

if param.scatter
    hold on
    plot(x,y,'.','color','k');
    hold off
end

end