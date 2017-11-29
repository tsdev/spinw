function varargout = activefigure(varargin)
% returns the handle of the active swplot figure
% 
% ### Syntax
% 
% `hFigure = swplot.activefigure`
%
% `swplot.activefigure(hFigure)`
% 
% ### Description
% 
% `hfigure = swplot.activefigure` returns the handle of the active swplot
% figure and makes it selected. If no swplot figure exists, the function
% throws an error.
%  
% `swplot.activefigure(hFigure)` makes the figure of `hFigure` handle the
% active figure.
% 
% ### Input Arguments
% 
% `hFigure`
% : Figure handle or figure number.
% 
% ### Output Arguments
% 
% `hFigure`
% : Handle of the active swplot figure.
% 
% ### See Also
% 
% [swplot.figure]
%

hFigure  = [];

% for plots, when there is no active figure, new figure will be created.
% For settings (such as legend, etc.) when there is no active figure, the
% last inactive figure will be used.
plotmode = false;

switch nargin
    case 1
        if ischar(varargin{1})
            plotmode = strcmpi(varargin{1},'plot');
        else
            hFigure = varargin{1};
        end
    case 2
        if ischar(varargin{1})
            varargin = varargin([2 1]);
        end
        hFigure = varargin{1};
        plotmode = strcmpi(varargin{2},'plot');
end

% tag for active figure
pref = swpref;
activeTag = pref.tag;

% tag for inactive figures
inactiveTag = ['inactive_' activeTag];

if isempty(hFigure)
    % find active figure
    hFigure = findobj('tag',activeTag);
    
    if isempty(hFigure)
        % no active figure create new figure for plotting to avoid overplot
        % on inactive figures
        if plotmode
            hFigure = swplot.figure;
        else
            % find the first inactive figure
            hFigure = findobj('tag',inactiveTag);
        end
        
        if isempty(hFigure)
            error('activefigure:NoFig','There is no swplot figure, use swplot.figure() to create a new window!')
            %hFigure = swplot.figure;
        end
        
        swplot.activefigure(hFigure(1));
        %warning('activefigure:Ativate','There is no active figure, activating the last used one!');
    end
    hFigure = hFigure(1);
    
    % make figure live
    %figure(hFigure);
else
    
    if ~(isnumeric(hFigure) || isgraphics(hFigure))
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