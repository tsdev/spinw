function hFigureOut = figure(mode)
% creates swplot figure
% 
% ### Syntax
% 
% `hFigure = swplot.figure`
%
% `hFigure = swplot.figure(mode)`
% 
% ### Description
% 
% `hFigure = swplot.figure` creates an empty figure with all the controls
% for modifying the plot and the 3D roation engine initialized that rotates
% the objects on the figure instead of the viewport. To plot anything onto
% the figure, the handle of the graphics object (after creating it using
% [matlab.surf], [matlab.patch], etc.) has to be added to the figure using
% the function [swplot.add].
%
% `hFigure = swplot.figure(mode)` defines settings for the figure.
%
% ### Input Arguments
% 
% `mode`
% : Optional string. If `'nohg'`, then no hgtransform object will be
%   used for fine object rotation. Can be usefull for certain
%   export functions, that are incompatible with [matlab.hgtransform]
%   objects. Default value is `'hg'` to use hgtransform.
% 
% ### See Also
% 
% [swplot.add] \| [matlab.hgtransform]
%

% how to avoid overplotting?
% if ~isempty(get(0,'CurrentFigure'))
%     oldAxis = get(gcf,'CurrentAxes');
% end
pref = swpref;

if nargin == 0
    mode = 'hg';
end

% Create new figure.
hFigure = figure;
% show hidden child handles, such as hToolbar
showHidden = get(0,'Showhidden');
set(0,'Showhidden','on')
set(gcf,'color','w')

% create drawing axis
hAxis = axes('Parent',hFigure);

% set camera to full manual control
set(hAxis,...
    'Position',         [0 0 1 1],...
    'Color',            'none',...
    'Box',              'off',...
    'Clipping',         'off',...
    'CameraPosition',   [0 0 1e4],...
    'CameraTarget',     [0 0 0],...
    'CameraUpVector',   [0 1 0],...
    'CameraViewAngle',  0.6,...
    'NextPlot',         'add',...
    'Tag',              'swaxis',...
    'XLim',             [-1000 1000],...
    'YLim',             [-1000 1000],...
    'ZLim',             [-1000 1000],...
    'Visible','off');

% fix 1:1 ratio for all scale
daspect([1 1 1]);
pbaspect([1 1 1]);

set(hFigure,...
	'Renderer',     'opengl',...
    'Name',         'swplot',...
    'DockControls', 'off',...
    'PaperType',    'A4',...
    'Tag',          pref.tag,...
    'Toolbar',      'figure',...
    'DeleteFcn',    @closeSubFigs...
    );

hToolbarT = get(hFigure,'children');

if verLessThan('matlab','8.4.0')
    hToolbar = hToolbarT(end);
else
    idx = 1;
    while ~isa(hToolbarT(idx),'matlab.ui.container.Toolbar') && idx<numel(hToolbarT)
        idx = idx + 1;
    end
    if idx>numel(hToolbarT)
        hToolbar = hToolbarT(end);
    else
        hToolbar = hToolbarT(idx);
    end
end

hButton  = get(hToolbar,'children');
switch mode
    case 'hg'
        delete(hButton([1:(end-3) end-1 end]));
    case 'nohg'
        % keep the object rotation button
        delete(hButton([1:(end-9) (end-7):(end-3) end-1 end]));
end

% Loads the icon variable, that contains all the used icons.
iconPath = [sw_rootdir 'dat_files' filesep 'icons.mat'];
try
    load(iconPath,'icon');
catch err
    error('figure:FileNotFound',['Icon definition file not found: '...
        regexprep(iconPath,'\' , '\\\') '!']);
end

% Change button backgrounds to the toolbar background color using
% undocumented features
try
    drawnow;
    hToolbar = findall(gcf,'tag','FigureToolBar');
    jToolbar = get(get(hToolbar,'JavaContainer'),'ComponentPeer');
    jColor   = jToolbar.getParent.getParent.getBackground();
    bkgColor = [get(jColor,'Red') get(jColor,'Green') get(jColor,'Blue')];
    fNames   = fieldnames(icon);
    
    for ii = 1:length(fNames)
        matSel = icon.(fNames{ii})*255;
        matSelR = matSel(:,:,1);
        matSelG = matSel(:,:,2);
        matSelB = matSel(:,:,3);
        idx = (matSelR == 212) & (matSelG == 208) & (matSelB == 200);
        matSelR(idx) = bkgColor(1);
        matSelG(idx) = bkgColor(2);
        matSelB(idx) = bkgColor(3);
        matSel = cat(3,matSelR,matSelB,matSelG);
        icon.(fNames{ii}) = matSel/255;
    end
catch %#ok<CTCH>
    warning('figure:JavaSupport','Some functionality won''t work, due to Java incompatibility');
end

if strcmp(mode,'hg')
    % only works in hg mode
    button.aAxis = uipushtool(hToolbar,'CData',icon.aAxis,'TooltipString','View direction: a axis',...
        'ClickedCallback',@(i,j)swplot.view('a',hFigure),'Separator','on');
    button.bAxis = uipushtool(hToolbar,'CData',icon.bAxis,'TooltipString','View direction: b axis',...
        'ClickedCallback',@(i,j)swplot.view('b',hFigure));
    button.cAxis = uipushtool(hToolbar,'CData',icon.cAxis,'TooltipString','View direction: c axis',...
        'ClickedCallback',@(i,j)swplot.view('c',hFigure));
    button.anyAxis = uipushtool(hToolbar,'CData',icon.anyAxis,'TooltipString','Set view direction',...
        'ClickedCallback','');
end

button.zoomOut = uipushtool(hToolbar,'CData',icon.zoomOut,'TooltipString','Zoom out',...
    'ClickedCallback',@(i,j)swplot.zoom(1/sqrt(2),hFigure),'Separator','on');
button.zoomIn  = uipushtool(hToolbar,'CData',icon.zoomIn,'TooltipString','Zoom in ',...
    'ClickedCallback',@(i,j)swplot.zoom(sqrt(2),hFigure));
button.setRange = uipushtool(hToolbar,'CData',icon.setRange,'TooltipString','Set range',...
    'ClickedCallback',{@swplot.setrangegui hFigure},'Separator','on');
button.figActive = uipushtool(hToolbar,'CData',icon.active,'TooltipString','Keep figure',...
    'Separator','on');
button.ver = uipushtool(hToolbar,'CData',icon.ver,'TooltipString','Show SpinW version',...
    'Separator','on','ClickedCallback',{@swplot.logo ''});
set(button.figActive,'ClickedCallback',{@activatefigure 'toggle'});

% activate figure
activatefigure([], [], 'initialize')

if strcmp(mode,'hg')
    hRotate = hgtransform('Parent',hAxis);
    hTranslate = hgtransform('Parent',hRotate);
    setappdata(hFigure,'h',hTranslate);
else
    hTranslate = gobjects(0);
    setappdata(hFigure,'h',hTranslate);
end

% save data to figure
setappdata(hFigure,'button',button);
setappdata(hFigure,'axis',hAxis);
setappdata(hFigure,'legend',struct('handle',gobjects(0),'text',{''},'type',[],'color',[],'name',{''}));
% add light
hLight = camlight('right');
set(hLight,'Parent',hAxis);
setappdata(hFigure,'light',hLight);
setappdata(hFigure,'icon',icon);
setappdata(hFigure,'base',eye(3));
swplot.zoom('auto',hFigure);
setappdata(hFigure,'tooltip',struct('handle',gobjects(0)));
%swplot.tooltip('off',hFigure);

% create empty object to store plot data
fNames = {'handle' 'number' 'name' 'type' 'label' 'position' 'text' 'legend' 'data'};
c0 = cell(1,0);
sInit = [fNames; repmat({c0},[1 numel(fNames)])];
setappdata(hFigure,'objects',struct(sInit{:}));

if nargout > 0
    hFigureOut = hFigure;
end

% restore showHidden property
set(0,'Showhidden',showHidden);

% add mouse control
swplot.mouse(hFigure);

% show figure
set(hFigure,'Visible','on');


    function activatefigure(~,~,mode)
        % activate/deactivate figure
        %
        % mode = initialize/toggle
        %
        
        inact = 'inactive_';
        
        fTag = get(hFigure,'tag');
        
        switch mode
            case 'initialize'
                fTag = [inact fTag];
            case 'toggle'
        end
        
        
        if ~isempty(strfind(fTag,inact)) || strcmp(mode,'initialize')
            % activate figure
            newTag = fTag((numel(inact)+1):end);
            set(hFigure,'tag',newTag);
            set(button.figActive,'CData',icon.active);
            hList = findobj('Tag',newTag);
            hList(hList==hFigure) = [];
            % The handle of the activate button is Figure.Child(end).Child(1).
            for jj = 1:numel(hList)
                bt = getappdata(hList(jj),'button');
                set(bt.figActive,'CData',icon.inactive);
                set(hList(jj),'Tag',[inact newTag]);
            end
            
        else
            % inactivate figure
            set(hFigure,'tag',[inact fTag]);
            set(button.figActive,'CData',icon.inactive);
        end
        
        
    end

    function closeSubFigs(~, ~)
        % close all sub figures on plot window closing
        
        % close the "Set Range" window if it is open
        if verLessThan('matlab','8.4.0')
            sub{1} = findobj('Tag',['setRange_' num2str(hFigure)]);
            sub{2} = findobj('Tag',['tooltipfig_' num2str(hFigure)]);
        else
            sub{1} = findobj('Tag',['setRange_' num2str(hFigure.Number)]);
            sub{2} = findobj('Tag',['tooltipfig_' num2str(hFigure.Number)]);
        end
        
        sub = sub(cellfun(@(C)~isempty(C),sub));
        close([sub{:}]);
        
        % delete figure
        delete(hFigure);
        
    end
end