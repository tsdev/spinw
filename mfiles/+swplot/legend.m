function legend(switch0,hFigure)
% draws legend to the swplot figure
%
% SWPLOT.LEGEND({switch}, {hFigure})
%
% Input:
%
% switch        String, can be 'on' or 'off' to switch the legend on/off.
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

switch switch0
    case 'on'
        switchon = true;
    case 'off'
        switchon = false;
    otherwise
        error('legend:WrongInput','Use on/off to switch legend!')
end

if ~switchon && ~isempty(lDat.handle) 
    % remove legend
    delete(lDat.handle)
    lDat.handle = gobjects(0);
    setappdata(hFigure,'legend',lDat);
    return
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
lAxis = axes(hFigure,'Units','pixel','Position',[dA dA lWidth+dA lHeight+dA],...
    'Visible','off','XLim',[0 lWidth],'YLim',[0 lHeight],'NextPlot','add');
lDat.handle = lAxis;

handle = rectangle('Position',[1 1 lWidth-2 lHeight-2],'FaceColor','w');

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
            handle(end+1) = rectangle('Position',[5 lHeight-20*ii+6 sRadius*2 sRadius],...
                'FaceColor',lColor(:,ii),'EdgeColor','k'); %#ok<*AGROW>
            %handle(end+1) = text(30,(lHeight-20*ii+10),lText{ii},'fontSize',fontSize,'color','k');
        case 2
            % dashed rectangle
            handle(end+1) = rectangle('Position',[6 lHeight-20*ii+6.9 sRadius*2/3-1 sRadius-1.5],...
                'FaceColor',lColor(:,ii),'EdgeColor',lColor(:,ii));
            handle(end+1) = rectangle('Position',[5+sRadius*2/3*2 lHeight-20*ii+6.9 sRadius*2/3-1 sRadius-1.5],...
                'FaceColor',lColor(:,ii),'EdgeColor',lColor(:,ii));
            handle(end+1) = rectangle('Position',[5 lHeight-20*ii+6 sRadius*2 sRadius],...
                'FaceColor','none','EdgeColor','k');
            %handle(end+1) = text(30,(lHeight-20*ii+10),lText{ii},...
            %    'fontSize',fontSize,'color','k');
        case 3
            % sphere
            handle(end+1) = swplot.ellipsoid([5+sRadius lHeight-20*ii+3+sRadius 0],eye(3)*sRadius,3);
    end
    % add text
    handle(end+1) = text(30,(lHeight-20*ii+10),lText{ii},'fontSize',fontSize,'color','k');
end

if any(lType==3)
    % add light
    camlight('right');
end

% take care that one cannot activate axis by clicking on it
set(handle,'Tag','legend');
set(handle,'hittest','off');
set(lAxis,'hittest','off','NextPlot','replace','Visible','on')

% save new handle
setappdata(hFigure,'legend',lDat);

% switch do plot axis where all 3D objects are
axes(getappdata(hFigure,'axis'));

end