function tooltip(text0,hFigure)
% creates tooltip axis on swplot figure
%
% SWPLOT.TOOLTIP({text}, {hFigure})
%
% Input:
%
% text          String, if it is 'on'/'off' the tooltip will be switched
%               on/off. Otherwise the text will be shown in the tooltip.
%               Default is 'on'. 
% hFigure       Handle of the swplot figure. Default is the selected
%               figure.
%
% Example:
%   swplot.figure
%   swplot.addcircle([0 0 0],[0 0 1],1)
%   swplot.tooltip
%

if nargin == 0
    text0 = 'on';
end

fontSize = swpref.getpref('fontsize',[]);

if nargin < 2
    % find active figure
    hFigure = swplot.activefigure;
end

tDat = getappdata(hFigure,'tooltip');

if isempty(tDat.handle)
    % create tooltip axis
    tDat.handle    = axes(hFigure,'Units','Normalized','Position',[0.01 0.9 0.1 0.2],'Visible','off');
    tDat.handle(2) = text(0.05,0.45,'','units','normalized','horizontalalignment','left',...
        'FontSize',fontSize,'VerticalAlignment','top');
    % avoid object get active for strange zooming effect
    set(tDat.handle,'hittest','off','tag','tooltip');
    
    % switch to plot axis where all 3D objects are
    axes(getappdata(hFigure,'axis'));

    % save new handles
    setappdata(hFigure,'tooltip',tDat);

end

switch text0
    case 'on'
        set(tDat.handle(2),'Visible','on');
    case 'off'
        set(tDat.handle(2),'Visible','off');
    otherwise
        set(tDat.handle(2),'String',text0);
end

end