function varargout = plot(obj, varargin)
% plots 3D model
% 
% ### Syntax
% 
% `plot(obj,Name,Value)`
% `hFigure = plot(obj,Name,Value)`
% 
% ### Description
% 
% `plot(obj,Name,Value)` plots the atoms and couplings stored in `obj` onto
% an [swplot] figure (see [swplot.figure]). The generated 3D plot can be
% rotated using the mouse and panning works by keeping the Ctrl/Control
% button pressed. There is information about every object on the figure
% (here called tooltips) that is shown when clicked on the object. The 3D
% view direction can be changed programatically using [swplot.view] while
% translations are controlled using the [swplot.translate]. Arbitrary
% transformation (combination of rotation and translation) can be
% introduced using the [swplot.transform]. All these transformation act as
% a global transformation, relative transformation of the 3D objects is
% only possible at creation by defining the transformed coordinates.
%  
% The `spinw.plot` function calls several high level plot routines to draw
% the different types of objects: [swplot.plotatom] (atoms),
% [swplot.plotmag] (magnetic moments), [swplot.plotion] (single ion
% properties), [swplot.plotbond] (bonds), [swplot.plotbase] (basis vectors)
% and [swplot.plotcell] (unit cells).
%  
% The high level `spinw.plot` function can send send parameters to any of
% the above plot group functions. The paramer name has to be of the format:
% `['plot group name' 'group option']`. For example to set the `color` option
% of the cell (change the color of the unit cell) use the option
% 'cellColor'. In this case `spinw.plot` will call the [swplot.plotcell]
% function with the `color` parameter set to the given value. For all the
% possible group plot function options see the corresponding help.
%  
% It is possible to switch off calling any of the subfunctions by using the
% option `['plot group name' 'mode']` set to `'none'`. For example to skip
% plotting of the atoms set the `'atomMode'` parameter to `'none'`:
% `spinw.plot('atomMode','none')`.
%  
% ### Name-Value Pair Arguments
%  
% These are global options, that each plot group function recognizes, these global
% options can be added without the group name.
% 
% `'range'`
% : Plotting range of the lattice parameters in lattice units,
%   in a matrix with dimensions of $[3\times 2]$. For example to plot the
%   first unit cell, use: `[0 1;0 1;0 1]`. Also the number unit cells can
%   be given along the $a$, $b$ and $c$ directions, e.g. `[2 1 2]`, this is
%   equivalent to `[0 2;0 1;0 2]`. Default value is the single unit cell.
% 
% `'unit'`
% : Unit in which the range is defined. It can be the following
%   string:
%   * `'lu'`        Lattice units (default).
%   * `'xyz'`       Cartesian coordinate system in \\ang units.
% 
% `'figure'`
% : Handle of the [swplot] figure. Default is the active figure.
% 
% `'legend'`
% : Whether to add legend to the plot, default value is `true`, for details
%   see [swplot.legend].
% 
% `'fontSize'`
% : Font size of the atom labels in pt units, default value is stored in
%   `swpref.getpref('fontsize')`.
% 
% `'nMesh'`
% : Resolution of the ellipse surface mesh. Integer number that is
%   used to generate an icosahedron mesh with $n_{mesh}$ number of
%   additional triangulation, default value is stored in
%   `swpref.getpref('nmesh')`.
% 
% `'nPatch'`
% : Number of points on the curve for the arrows and cylinders,
%   default value is stored in `swpref.getpref('npatch')`.
% 
% `'tooltip'`
% : If `true`, the tooltips will be shown when clicking on the plot.
%   Default value is `true`.
% 
% `'shift'`
% : Column vector with 3 elements, all objects will be shifted by
%   the given value. Default value is `[0;0;0]`.
% 
% `'replace'`
% : Replace previous plot if `true`. Default value is `true`.
% 
% `'translate'`
% : If `true`, all plot objects will be translated to the figure
%   center. Default is `true`.
% 
% `'zoom'`
% : If `true`, figure will be automatically zoomed to the ideal scale.
%   Default value is `true`.
%
% ### See Also
%  
% [swplot.plotatom] \| [swplot.plotmag] \| [swplot.plotion] \| 
% [swplot.plotbond] \| [swplot.plotbase] \| [swplot.plotcell]
%
pref = swpref;
fid = pref.fid;

% preparation
fprintf0(fid,'Creating 3D plot... \n');

if numel(varargin) == 1
    % handle input structures
    if ~isstruct(varargin{1})
        error('spinw:plot:WrongOptions','Expected option-value pairs')
    end
    % convert struct to cell
    input = [fieldnames(varargin{1}) struct2cell(varargin{1})]';
    input = input(:)';
else
    input = varargin;
end

% check varargin
if mod(numel(input),2) ~= 0
    error('spinw:plot:WrongOptions','Expected option-value pairs')
end

% convert input into a struct
optName = lower(input(1:2:end));
optVal  = input(2:2:end);

% find global options that won't be called
inpForm.fname  = {'tooltip' 'figure' 'zoom' 'translate' 'replace'};
inpForm.defval = {true      []       true   true        true     };
inpForm.size   = {[1 1]     [1 1]    [1 1]  [1 1]       [1 1]    };
inpForm.soft   = {false     true     false  false       false    };

noName = inpForm.fname;
selNoG = ismember(optName,noName);
optNameNoG = optName(selNoG);
optValNoG  = optVal(selNoG);

optName = optName(~selNoG);
optVal  = optVal(~selNoG);

optNo = [optNameNoG;optValNoG];
param = sw_readparam(inpForm, optNo{:});

if isempty(param.figure)
    hFigure = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

% find global options:
globalName = {'range' 'unit' 'legend' 'fontsize' 'nmesh' 'npatch' 'shift' 'copy'};
selG = ismember(optName,globalName);
% take care of global figure handle
optNameG = [optName(selG) {'figure' 'tooltip' 'zoom' 'translate' 'replace'}];
optValG  = [optVal(selG)  {hFigure  false     false   false      false    }];
% all global option
varargG  = [optNameG;optValG];

optName = optName(~selG);
optVal  = optVal(~selG);

% check some global option for intervention
range0 = [0 1;0 1;0 1];

inpForm.fname  = {'range' 'unit'};
inpForm.defval = {range0  'lu'  };
inpForm.size   = {[-1 -2] [1 -3]};
inpForm.soft   = {false   false };

warn0 = warning('query');
warning('off','sw_readparam:UnreadInput')
warning('off','plotbond:EmptyPlot')

paramG = sw_readparam(inpForm, varargG{:});

% sort options according to plot name
plotName = {'atom' 'mag' 'bond' 'ion' 'cell' 'base' 'chem'};
nFun     = numel(plotName);
plotIdx  = false(nFun,numel(optName));
optShort = optName;

if ~isempty(optName)
    for ii = 1:nFun
        fIdx = strfind(optName,plotName{ii});
        plotIdx(ii,:) = cellfun(@(C)~isempty(C) && C(1)==1,fIdx);
        optShort(plotIdx(ii,:)) = cellfun(@(C)C((numel(plotName{ii})+1):end),optName(plotIdx(ii,:)),'UniformOutput',false);
    end
end

if ~all(any(plotIdx,1))
    errIdx = find(~any(plotIdx,1));
    error('spinw:plot:WrongInput',['Wrong option: ''' optName{errIdx(1)} '''!']);
end

% list of plot functions
plotFun = {@swplot.plotatom @swplot.plotmag @swplot.plotbond @swplot.plotion @swplot.plotcell @swplot.plotbase @swplot.plotchem};

% check switches (mode == 'none'), default is to show plot
switchFun = true(1,nFun);
% default is not to call swplot.plotchem()
switchFun(end) = false;

for ii = 1:nFun
    if isempty(optShort) || isempty(plotIdx)
        modeIdx = [];
    else
        modeIdx = find(ismember(optShort,'mode') & plotIdx(ii,:));
    end
    
    if numel(modeIdx)>1
        error('spinw:plot:WrongInput',['Muiltiple options: ''' optName{modeIdx(1)} '''!']);
    end
    
    if numel(modeIdx) == 1
        switchFun(ii) = ~strcmp(optVal{modeIdx},'none');
    end
    
    switch plotName{ii}
        % additional tests
        case 'atom'
        case 'mag'
            if isempty(obj.mag_str.F)
                % don't plot magnetic structure
                switchFun(ii) = false;
            end
        case 'bond'
            if isempty(obj.coupling.atom1)
                % don't plot bonds
                switchFun(ii) = false;
            end
        case 'ion'
        case 'cell'
        case 'base'
    end
    
end

if ~any(switchFun)
    % plot nothing
    return
end

% only call selected functions
plotFun = plotFun(switchFun);
plotIdx = plotIdx(switchFun,:);

% short arguments
vararg  = [optShort;optVal];

warning('on','sw_readparam:UnreadInput')

% remove objects with the same name if necessary
if param.replace
    % delete all object that is created by spinw.plot
    % find objects to be deleted
    name0 = {'base' 'base_label' 'atom' 'atom_label' 'bond' 'bond_mat' 'cell' 'mag' 'ion' 'ion_edge' 'chem'};
    sObj = swplot.findobj(hFigure,'name',name0);
    % delete them!
    swplot.delete([sObj(:).number]);
    
    % remove old legend entries
    lDat = getappdata(hFigure,'legend');
    if ~isempty(lDat.name)
        lIdx = ~ismember(lDat.name,{'atom' 'bond'});
        lDat.color = lDat.color(:,lIdx);
        lDat.type  = lDat.type(:,lIdx);
        lDat.name  = lDat.name(:,lIdx);
        lDat.text  = lDat.text(:,lIdx);
        setappdata(hFigure,'legend',lDat);
        % redraw legend
        swplot.legend('refresh',hFigure);
    end
end

% call plot functions
for ii = 1:numel(plotFun)
    selArg = vararg(:,plotIdx(ii,:));
    % add any global parameter that might be missing
    if ~isempty(selArg)
        selG = ~ismember(optNameG,selArg(1,:));
        selArg = [varargG(:,selG) selArg]; %#ok<AGROW>
    else
        selArg = varargG;
    end
    
    if ii == 1
        selArg = [selArg {'obj'; obj}]; %#ok<AGROW>
    end
    
    hFigure = plotFun{ii}(selArg{:});
    
end

if param.tooltip
    swplot.tooltip('on',hFigure);
end

% save all options
setappdata(hFigure,'param',input);

% change range, if the number of unit cells are given
if numel(paramG.range) == 3
    paramG.range = [ zeros(3,1) paramG.range(:)];
elseif numel(paramG.range) ~=6
    error('spinw:plot:WrongInput','The given plotting range is invalid!');
end

% save plotting range and unit
rDat.range = paramG.range;
rDat.unit  = paramG.unit;
setappdata(hFigure,'range',rDat);

if nargout > 0
    varargout{1} = hFigure;
end

if param.zoom
    swplot.zoom('auto',hFigure);
end

if param.translate
    swplot.translate('auto',hFigure);
end

if param.tooltip
    swplot.tooltip('on',hFigure);
end

set(hFigure,'Name', 'SpinW plot()');
drawnow

% patch handles
hPatch = findobj('type','patch');
nFaces    = sum(cellfun(@(C)size(C,1),get(hPatch,'Faces')));
nVertices = sum(cellfun(@(C)size(C,1),get(hPatch,'Vertices')));

% ready
fprintf0(fid,'...%dk faces and %dk vertices are drawn!\n',round(nFaces/1e3),round(nVertices/1e3));

warning(warn0);

end