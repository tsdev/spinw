function varargout = legend(switch0,hFigure)
% draws legend to the swplot figure
%
% SWPLOT.LEGEND({switch}, {hFigure})
%
% status = SWPLOT.LEGEND
%
% Input:
%
% switch        One of the following string:
%                   'on'                show legend,
%                   'off'               hide legend,
%                   'refresh'           redraw legend,
%                   {'-','--','none'}   change the linestyle of the legend
%                                       frame.
%               Default is 'on'.
% hFigure       Handle of the swplot figure. Default is the selected
%               figure.
%
% Example:
%   swplot.figure
%   swplot.addcircle([0 0 0],[0 0 1],1)
%   swplot.legend
%

if nargin == 0
    switch0 = 'on';
end

if nargin < 2
    % find active figure
    hFigure = swplot.activefigure;
end

lDat = getappdata(hFigure,'legend');

if nargout > 0
    if isempty(lDat.handle)
        varargout{1} = 'off';
    else
        varargout{1} = 'on';
    end
    return
end

refresh = false;

switch switch0
    case 'on'
        switchon = true;
    case 'off'
        switchon = false;
    case 'refresh'
        switchon = ~isempty(lDat.handle);
        refresh  = true;
    case {'-' '--' 'none'} %'frame'
        switchon = true;
    otherwise
        error('legend:WrongInput','Use on/off to switch legend!')
end

if (~switchon && ~isempty(lDat.handle)) || refresh
    % remove legend
    delete(lDat.handle)
    lDat.handle = gobjects(0);
    setappdata(hFigure,'legend',lDat);
    if ~refresh
        return
    end
end

if isempty(lDat.type)
    % no legend entries don't show anything
    return
end

% dimensions of legend
lHeight = numel(lDat.type)*20;
lWidth  = 80;
sRadius = 8;
dA      = 2;

if ~isempty(lDat.handle)
    delete(lDat.handle)
end

% create new axis for legend
hAxis = axes('Parent',hFigure,'Units','pixel','Position',[dA dA lWidth+dA lHeight+dA],...
    'Visible','off','XLim',[0 lWidth],'YLim',[0 lHeight+2],'NextPlot','add','box','off');

hObject = rectangle('Parent',hAxis,'Position',[1 1 lWidth-2 lHeight],'FaceColor','w');

lDat.handle = [hAxis hObject];

% extract legend data
lType  = lDat.type;
lText  = lDat.text;
lColor = lDat.color;

% get the stored fontsize
fontSize = swpref.getpref('fontsize',[]);

for ii = 1:numel(lType)
    switch lType(ii)
        case 1
            % colored rectangle
            hObject(end+1) = rectangle('Parent',hAxis,'Position',[5 lHeight-20*ii+6 sRadius*2 sRadius],...
                'FaceColor',lColor(:,ii),'EdgeColor','k'); %#ok<*AGROW>
        case 2
            % dashed rectangle
            hObject(end+1) = rectangle('Parent',hAxis,'Position',[6 lHeight-20*ii+6.9 sRadius*2/3-1 sRadius-1.5],...
                'FaceColor',lColor(:,ii),'EdgeColor',lColor(:,ii));
            hObject(end+1) = rectangle('Parent',hAxis,'Position',[5+sRadius*2/3*2 lHeight-20*ii+6.9 sRadius*2/3-1 sRadius-1.5],...
                'FaceColor',lColor(:,ii),'EdgeColor',lColor(:,ii));
            hObject(end+1) = rectangle('Parent',hAxis,'Position',[5 lHeight-20*ii+6 sRadius*2 sRadius],...
                'FaceColor','none','EdgeColor','k');
        case 3
            % sphere
            hObject(end+1) = swplot.ellipsoid(hAxis,[5+sRadius lHeight-20*ii+3+sRadius 0],eye(3)*sRadius,3);
            set(hObject(end),'FaceColor',lColor(:,ii));
    end
    % add text
    hObject(end+1) = text(30,(lHeight-20*ii+10),lText{ii},'Parent',hAxis,'fontSize',fontSize,'color','k');
end

hText = hObject(ismember(get(hObject,'Type'),'text'));
% size of text
tSize = get(hText,'Extent');
if iscell(tSize)
    tSize = reshape([tSize{:}],4,[])';
end
tSize = max(sum(tSize(:,[1 3]),2))+5;
% extend rectangle and axis
set(hObject(1),'Position',[1 1 tSize lHeight]);
set(hAxis,'Position',[dA dA tSize+dA lHeight+dA],'XLim',[0 tSize+1]);

if any(lType==3)
    % add light
    hLight = camlight('right');
    set(hLight,'Parent',hAxis);
    material(hAxis,'shiny');
end

if ismember(switch0,{'-' '--' 'none'}) && ~isempty(lDat.handle)
    % change rectangle line style
    set(lDat.handle(2),'LineStyle',switch0);
end

% take care that one cannot activate axis by clicking on it
set(hObject,'Tag','legend');
set(hObject,'hittest','off');
set(hAxis,'hittest','off','NextPlot','replace','Visible','off')

% save new handle
setappdata(hFigure,'legend',lDat);

% switch to plot axis where all 3D objects are
%axes(getappdata(hFigure,'axis'));

end