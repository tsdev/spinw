function udat = spaghetti(dat,varargin)
% creates spaghetti plot from d2d objects
%
% udat = NDBASE.SPAGHETTI(dat, 'option1', value1, ...)
%
% Input:
%
% dat       Array of d2d class objects with nDat member. Typically Horace
%           cuts.
%
% Options:
%
% flip      Vector of nDat logical values, if any value is true the
%           corresponding data is flipped along the horizontal axis.
%           Default is false(1,nDat).
% label     Cell of strings for the x-axis label. There has to be nDat+1
%           strings.
% dashed    If true, black dashed line is drawn between cuts. Default is
%           true.
% ylim      If true, the upper limit of the vertical axis is automatically
%           determined. Default is true.
% plot      If true, plot is created. Default is true.
% pad       If true, the output data is padded with nans on the top and
%           right side. Then it can be plotted directly using surf():
%               udat = ndbase.spaghetti(dat,'pad',true);
%               coo = cell(1,2);
%               [coo{:}] = ndgrid(udat.x,udat.y);
%               surf(coo{:},udat.sig);
%
% Output:
%
% udat      United data in simple struct with the following fields:
%               x       Column vector of x-coordinate edge bins (nX+1).
%               y       Column vector of y-coordinate edge bins (nY+1).
%               sig     Signal matrix with dimensions [nX nY].
%               err     Error matrix with dimensions [nX nY].
%               ylabel  Label of the y-axis, string.
%               xlim    Limits of the x-axis, row vector with 2 elements.
%               ylim    Limits of the y-axis, row vector with 2 elements.
%               clim    Limits of the c-axis, row vector with 2 elements.
%

inpForm.fname  = {'flip' 'label' 'dashed' 'ylim' 'plot' 'pad'};
inpForm.defval = {[]     []      true     'auto' true   false};
inpForm.size   = {[1 -1] [1 -2]  [1 1]    [1 -3] [1 1]  [1 1]};
inpForm.soft   = {true   true    false    false  false  false};

param = sw_readparam(inpForm, varargin{:});

% number of data pieces to plot
nDat = numel(dat);

if isempty(param.flip)
    param.flip = false(1,nDat);
else
    if numel(param.flip)~=nDat
        error('spaghetti:WrongOption','The ''flip'' option should containt true/false for each data piece!');
    end
    param.flip = logical(param.flip);
end

if ~isa(dat,'d2d')
    error('spaghetti:WrongDataType','The given data is not d2d type!')
end

if isempty(param.label)
else
    if numel(param.label)~=(nDat+1)
        error('spaghetti:WrongOption','The ''label'' options should contain a string label for each x edge!')
    end
end

spags = struct(dat(1));
sig   = spags.s;
err   = spags.e;
% remove signal from empty channels
sig(spags.npix==0) = nan;
err(spags.npix==0) = nan;

if param.flip(1)
    sig = flipud(sig);
    err = flipud(err);
end

ax0   = spags.p{1}(:)*spags.ulen(spags.pax(1));
%ax{1} = ax0-mean(ax0(1:2));
ax{1} = ax0;
ax{2} = spags.p{2}(:);
dash  = zeros(1,nDat-1);

for ii = 2:nDat
    spags = struct(dat(ii));
    
    sig0 = spags.s(2:end,:);
    err0 = spags.e(2:end,:);
    npix = spags.npix(2:end,:);
    
    sig0(npix==0) = nan;
    err0(npix==0) = nan;
    
    if param.flip(ii)
        sig0 = flipud(sig0);
        err0 = flipud(err0);
    end
    
    sig = [sig; sig0]; %#ok<AGROW>
    err = [err; err0]; %#ok<AGROW>
    
    dash(ii-1) = mean(ax{1}(end-1:end));
    ax0   = spags.p{1}(2:end)*spags.ulen(spags.pax(1));
    % fit the axes together
    ax0   = ax0-ax0(1)+ax{1}(end);
    ax(1) = {[ax{1}; ax0(2:end)]};
    
    % check ax(2)
    if numel(spags.p{2})~=numel(ax{2}) || any(abs(spags.p{2}(:)-ax{2})>10*eps)
        error('spaghetti:WrongInput','The data has unequal vertical axes!')
    end
end

% axis limits
axLim = [mean(ax{1}(1:2)) mean(ax{1}(end-1:end)) mean(ax{2}(1:2)) mean(ax{2}(end-1:end))];
if strcmp(param.ylim,'auto')
    % auto y-limit
    idx = (find(any(~isnan(sig),1),1,'last'));
    if ~isempty(idx)
        axLim(3:4) = [mean(ax{2}(1:2)) mean(ax{2}(idx+[0 1]))];
    end
end

% color limits
if min(sig(~isnan(sig)))>0
    cLim = [0 max(sig(:))];
else
    cLim = [min(sig(:)) max(sig(:))];
end

% draw the plot
if param.plot
    [xx, yy] = ndgrid(ax{:});
    % pad signal
    sigP = sig;
    sigP(end+1,:) = nan;
    sigP(:,end+1) = nan;
    hAxis = gca;
    surf(xx,yy,xx*0,sigP,'EdgeAlpha',0,'Parent',hAxis);
    hold on
    view(2)
    
    if param.dashed && nDat>1
        xD = [dash(:) dash(:) nan*dash(:)]';
        yD = repmat([min(ax{2}) max(ax{2}) nan],[size(xD,2) 1])';
        line(xD(:),yD(:),'LineStyle','--','Color','k');
    end
    
    if ~isempty(param.label)
        hAxis.XTick = [0 dash mean(ax{1}(end-1:end))];
        hAxis.XTickLabel = param.label;
    end
    
    axis(axLim);
    caxis(cLim);
    colorbar
    
    % figure cosmetics
    ylabel(spags.ulabel(spags.pax(2)));
    box on
    grid off
    colormap(cm_viridis(500))
end

udat.x      = ax{1};
udat.y      = ax{2};
udat.sig    = sig;
udat.err    = sqrt(err);
udat.ylabel = spags.ulabel{spags.pax(2)};
udat.xlim   = axLim(1:2);
udat.ylim   = axLim(3:4);
udat.clim   = cLim;

if param.pad
    udat.sig(end+1,:) = nan;
    udat.sig(:,end+1) = nan;
end

end