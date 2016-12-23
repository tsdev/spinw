function tooltip(switch0,hFigure)
% creates tooltip axis on swplot figure
%
% SWPLOT.TOOLTIP({switch}, {hFigure})
%
% Input:
%
% switch        String, can be 'on' or 'off' to switch the tooltip on/off.
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
    switch0 = 'on';
end

fontSize = swpref.getpref('fontsize',[]);

if nargin < 2
    % find active figure
    hFigure = swplot.activefigure;
end

tDat = getappdata(hFigure,'tooltip');

if ~isempty(tDat.handle)
    % remove old tooltip
    delete(tDat.handle);
    tDat.handle = [];
end

switch switch0
    case 'on'
        tDat.handle    = axes(hFigure,'Units','Normalized','Position',[0.01 0.9 0.1 0.2],'Visible','off');
        tDat.handle(2) = text(0.05,0.45,'','units','normalized','horizontalalignment','left',...
            'FontSize',fontSize,'VerticalAlignment','top');
        % avoid object get active for strange zooming effect
        set(tDat.handle,'hittest','off','tag','tooltip');
    case 'off'
    otherwise
        error('tooltip:WrongInput','Use on/off to switch tooltip!')
end

% save new handles
setappdata(hFigure,'tooltip',tDat);

% switch do plot axis where all 3D objects are
axes(getappdata(hFigure,'axis'));

end