function hFigure = activefigure(hFigure)
% returns the handle of the active swplot figure
%
% hFigure = SWPLOT.ACTIVEFIGURE
%
% It gives the handle of the active swplot figure and make it selected. If
% no swplot figure is open, throws an error.
%
% SWPLOT.ACTIVEFIGURE(hFigure)
%
% It makes hFigure the active figure.
%
% Input:
%
% hFigure       Figure handle or figure number.
%
% Output:
%
% hFigure       Figure handle.
%
% See also SWPLOT.FIGURE.
%

defTag = swpref.getpref('tag',[]);

if nargin == 0
    % find active figure
    hFigure = findobj('tag',defTag);
    if isempty(hFigure)
        error('activefigure:NoFig','There is no swplot figure, use swplot.figure() first, to create a new window!')
    end
    hFigure = hFigure(1);
    
    % make figure live
    figure(hFigure);
else
    
    if ~isgraphics(hFigure) || ~isnumeric(hFigure)
        error('activefigure:WrongInput','The given input is not a handle of an swplot figure!');
    elseif isnumeric(hFigure)
        hFigure = figure(hFigure);
    end
    
    % activate figure
    inact = 'inactive_';
    fTag = get(hFigure,'tag');
    
    if strcmp(fTag,defTag)
        % figure is already active, nothing to do
    elseif strcmp(fTag,[inact defTag])
        % figure is inactive
        % find active figure
        hFigureActive = findobj('tag',defTag);
        % deactivate
        bt = getappdata(hFigureActive,'button');
        icon = getappdata(hFigureActive,'icon');
        set(bt.figActive,'CData',icon.inactive);
        set(hFigureActive,'Tag',[inact defTag]);
        
        % activate new figure
        bt = getappdata(hFigure,'button');
        set(bt.figActive,'CData',icon.active);
        set(hFigure,'Tag',defTag);
        
    else
        error('activefigure:WrongInput','The given input is not a handle of an swplot figure!');
    end
end