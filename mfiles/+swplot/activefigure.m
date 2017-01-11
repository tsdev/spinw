function varargout = activefigure(mode,hFigure)
% returns the handle of the active swplot figure
%
% {hFigure} = SWPLOT.ACTIVEFIGURE
%
% It gives the handle of the active swplot figure and make it selected. If
% no swplot figure exists, throws an error.
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

if nargin == 0
    ishg = true;
    hFigure = [];
elseif ~ischar(mode)
    hFigure = mode;
    ishg    = true;
elseif nargin == 1
    hFigure = [];
    ishg = strcmp(mode,'hg');
elseif nargin == 2
    ishg = strcmp(mode,'hg');
end

if ishg
    mode = 'hg';
else
    mode = 'nohg';
end

% tag for active figure
activeTag = swpref.getpref('tag',[]);

% tag for inactive figures
inactiveTag = ['inactive_' activeTag];

if isempty(hFigure)
    % find active figure
    hFigure = findobj('tag',activeTag);
    % find the right type of figure
    hFigure = hFigure(swplot.ishg(hFigure)==ishg);
    
    if isempty(hFigure)
        % find the first inactive figure
        hFigure = findobj('tag',inactiveTag);
        hFigure = hFigure(swplot.ishg(hFigure)==ishg);
        
        if isempty(hFigure)
            %error('activefigure:NoFig','There is no swplot figure, use swplot.figure() to create a new window!')
            hFigure = swplot.figure(mode);
        end
        swplot.activefigure(hFigure(1));
        %warning('activefigure:Ativate','There is no active figure, activating the last used one!');
    end
    hFigure = hFigure(1);
    
    % make figure live
    %figure(hFigure);
else
    
    if ~(isgraphics(hFigure) || isnumeric(hFigure))
        error('activefigure:WrongInput','The given input is not a handle of an swplot figure!');
    elseif isnumeric(hFigure)
        hFigure = figure(hFigure);
    end
    
    % activate figure
    fTag = get(hFigure,'tag');
    
    if strcmp(fTag,activeTag)
        % figure is already active, nothing to do
    elseif strcmp(fTag,inactiveTag)
        % figure is inactive
        % find active figure
        hFigureActive = findobj('tag',activeTag);
        % deactivate
        if ~isempty(hFigureActive)
            bt = getappdata(hFigureActive,'button');
            icon = getappdata(hFigureActive,'icon');
            set(bt.figActive,'CData',icon.inactive);
            set(hFigureActive,'Tag',inactiveTag);
        end
        
        % activate new figure
        bt = getappdata(hFigure,'button');
        icon = getappdata(hFigure,'icon');
        set(bt.figActive,'CData',icon.active);
        set(hFigure,'Tag',activeTag);
        
    else
        error('activefigure:WrongInput','The given input is not a handle of an swplot figure!');
    end
    
end

if nargout>0
    varargout{1} = hFigure;
end

end