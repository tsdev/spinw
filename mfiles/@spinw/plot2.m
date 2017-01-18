function varargout = plot2(obj, varargin)
% plots crystal structure, magnetic structure, anisotropy and bonds
%
% PLOT(obj, 'Option1', Value1, ...)
%
% hFigure = PLOT(obj, 'Option1', Value1, ...)
%
% The function plots the atoms and couplings stored in obj onto an swplot
% figure.
%
% The function calls several plot routine to draw different group of
% objects to the figure: atom (atoms), mag (magnetic moments), ion (single
% ion properties), bond (bonds), base (basis vectors) and cell (unit
% cells).
%
% Each group has a corresponding plot function with its own input options,
% for the allowed options, please check the help of the functions below.
%
% atom  swplot.plotatom
% mag   swplot.plotmag
% ion   swplot.plotion
% bond  swpplot.plotbond
% base  swplot.plotbase
% cell  swplot.plotcell
%
% To provide an option to any of these sub functions add the name of the
% group to the option. For example to set the 'color' option of the cell
% (change the color of the unit cell) use the option 'cellColor'.
%
% It is possible to switch of any of the subfunctions by using the option
% [groupName 'mode'] set to 'none'. For example to skip plotting the atoms
% use the 'atomMode' option set to 'none'.
%
% There are also global options, that each group recognizes, these global
% options can be added without the group name.
%
% Global options:
%
% range     Plotting range of the lattice parameters in lattice units,
%           dimensions are [3 2]. For example to plot the first unit cell,
%           use: [0 1;0 1;0 1]. Also the number unit cells can be given
%           along the a, b and c directions: [2 1 2], that is equivalent to
%           [0 2;0 1;0 2]. Default is the single unit cell.
% rangeunit Unit in which the range is defined. It can be the following
%           string:
%               'lu'        Lattice units (default).
%               'xyz'       Cartesian coordinate system in Angstrom units.
% figure    Handle of the swplot figure. Default is the selected figure.
% legend    Whether to add the plot to the legend, default is true.
% fontSize  Font size of the atom labels in pt, default is stored in
%           swpref.getpref('fontsize').
% nMesh     Resolution of the ellipse surface mesh. Integer number that is
%           used to generate an icosahedron mesh with #mesh number of
%           additional triangulation, default value is stored in
%           swpref.getpref('nmesh')
% nPatch    Number of points on the curve for the arrows and cylinders,
%           default value is stored in swpref.getpref('npatch').
% tooltip   If true, the tooltips will be shown when clicking on the plot.
%           Default is true.
% shift     Column vector with 3 elements, all objects will be shifted by
%           the given value. Default value is [0;0;0].
% replace   Replace previous plots if true. Default is true.
%


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

% find global options that won't be calles
inpForm.fname  = {'tooltip' 'figure'};
inpForm.defval = {true      []      };
inpForm.size   = {[1 1]     [1 1]   };
inpForm.soft   = {false     true    };

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
globalName = {'range' 'rangeunit' 'legend' 'fontsize' 'nmesh' 'npatch' 'shift' 'replace'};
selG = ismember(optName,globalName);
% take care of global figure handle
optNameG = [optName(selG) {'figure'}];
optValG  = [optVal(selG) {hFigure}];
% all global option
varargG  = [optNameG;optValG];

optName = optName(~selG);
optVal  = optVal(~selG);

% check some global option for intervention
range0 = [0 1;0 1;0 1];

inpForm.fname  = {'range'};
inpForm.defval = {range0 };
inpForm.size   = {[-1 -2]};
inpForm.soft   = {false  };
warning('on','sw_readparam:UnreadInput')
warn0 = warning('query');
warning('off','sw_readparam:UnreadInput');
paramG = sw_readparam(inpForm, varargG{:});
warning(warn0);

% sort options according to plot name
plotName = {'atom' 'mag' 'bond' 'base' 'cell'}; % ion
nFun     = numel(plotName);
plotIdx  = false(nFun,numel(optName));
optShort = optName;

for ii = 1:nFun
    fIdx = strfind(optName,plotName{ii});
    plotIdx(ii,:) = cellfun(@(C)~isempty(C) && C(1)==1,fIdx);
    optShort(plotIdx(ii,:)) = cellfun(@(C)C((numel(plotName{ii})+1):end),optName(plotIdx(ii,:)),'UniformOutput',false);
end

if ~all(any(plotIdx,1))
    errIdx = find(~any(plotIdx,1));
    error('spinw:plot:WrongInput',['Wrong option: ''' optName{errIdx(1)} '''!']);
end

% check switches (mode == 'none')
switchFun = false(1,nFun);

for ii = 1:nFun
    modeIdx = find(ismember(optShort,'mode') & plotIdx(ii,:));
    if numel(modeIdx)>1
        error('spinw:plot:WrongInput',['Muiltiple options: ''' optName{modeIdx(1)} '''!']);
    end
    
    if numel(modeIdx) == 1
        switchFun(ii) = ~strcmp(optVal{modeIdx},'none');
    else
        switchFun(ii) = true;
    end
    
end

if ~any(switchFun)
    % plot nothing
    return
end

% only call selected functions
plotFun = {@swplot.plotatom @swplot.plotmag @swplot.plotbond @swplot.plotcell @swplot.plotbase}; % @swplot.plotion
plotFun = plotFun(switchFun);
plotIdx = plotIdx(switchFun,:);

% short arguments
vararg  = [optShort;optVal];

warning('on','sw_readparam:UnreadInput')

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
% save plotting range
setappdata(hFigure,'range',paramG.range);

if nargout > 0
    varargout{1} = hFigure;
end

end