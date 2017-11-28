function varargout = plotcell(varargin)
% plots unit cell
% 
% ### Syntax
% 
% `swplot.plotcell(Name,Value)`
% 
% `hFigure = swplot.plotcell(Name,Value)`
%
% ### Description
% 
% `swplot.plotcell(Name,Value)` plots the edges of the unit cell or
% multiple unit cells.
% 
% ### Name-Value Pair Arguments
% 
% `'mode'`
% : String that determines how the cells are plotted:
%   * `'inside'`    unit cells are plotted inside the volume defined by the
%                   given `range` parameter (default),
%   * `'single'`    a single unit cell is plotted at the origin,
%   * `'outside'`   unit cells are plotted inclusive the volume defined
%                   by the `range` parameter.
% 
% `'color'`
% : Color of the edges of the cells, one of the following values:
%   * `'auto'`      all edges will be black,
%   * `'colorname'` all edges will have the same color defined by the
%                   string, e.g. `'red'`,
%   * `[R G B]`     all edges will have the same color defined by the RGB
%                   code.
% 
% `'lineStyle'`
% : Determines the line style of the edges. Possible values are:
%   * `'--'`        dahsed edges (default),
%   * `'-'`         edges as continuous lines,
%   * `':'`         edges as dotted lines,
%   * `'-.'`        edges as dash-dotted lines.
% 
% `'lineWdith'`
% : Line width of the edges, default value is 1 pt.
% 
% #### General paraters
%
% These parameters have the same effect on any of the `swplot.plot...`
% functions.
%
% `'obj'`
% : [spinw] object.
% 
% `'unit'`
% : Unit in which the plotting range is defined. It can be one of the
%   following strings:
%   * `'lu'`        plot range is defined in lattice units (default),
%   * `'xyz'`       plot range is defined in the $xyz$ Cartesian coordinate
%                   system in \\ang units.
%
% `'range'`
% : Defines the plotting range. Depending on the `unit` parameter, the
%   given range can be in lattice units or in units of the $xyz$ Cartesian
%   coordinate system. It is either a matrix with dimensions of $[3\times
%   2]$ where the first and second columns define the lower and upper plot
%   limits respectively. It can be alternatively a vector with three
%   elements `[a,b,c]` which is equivalent to `[0 a;0 b;0 c]`. Default
%   value is `[0 1;0 1;0 1]` to show a single cell.
% 
% `'figure'`
% : Handle of the swplot figure. Default value is the active figure handle.
% 
% `'legend'`
% : Whether to show the plot on the legend, default value is `true`.
%
% `'tooltip'`
% : If `true`, the tooltips will be shown when clicking on the plot
%   objects. Default value is `true`.
% 
% `'shift'`
% : Column vector with 3 elements, all object positions will be
%   shifted by the given value in \\ang units. Default value is
%   `[0;0;0]`.
% 
% `'replace'`
% : If `true` the plot will replace the previous plot of the same type.
%   Default value is `true`.
% 
% `'translate'`
% : If `true`, the plot will be centered, independent of the range. Default
%   value is `false`.
% 
% `'zoom'`
% : If `true`, the swplot figure will be zoomed to make the plot objects
%   cover the full figure. Default is `true`.
% 
% `'copy'`
% : If `true`, a clone of the [spinw] object will be saved in the
%   swplot figure data which can be retwrived using
%   `swplot.getdata('obj')`. If `false`, the handle of the original [spinw]
%   object is saved which is linked to the input `obj` and so it changes
%   when `obj` is changed. Default value is `false`.
% 
% `nMesh`
% : Mesh of the ellipse surface, a triangulation class object or an
%   integer that used to generate an icosahedron mesh with $n_{mesh}$
%   number of additional subdivision into triangles. Default value is
%   stored in `swpref.getpref('nmesh')`, see also [swplot.icomesh].
% 
% `nPatch`
% : Number of vertices on any patch object that is not the icosahedron,
%   default value is stored in `swpref.getpref('npatch')`.
%
% ### Output Arguments
% 
% `hFigure`
% : Handle of the swplot figure.
%
% The name of the objects `'cell'`.
% To find the handles and the corresponding data on these objects, use e.g.
% sObject = swplot.findobj(hFigure,'name','cell')`.
%

% default values
nMesh0    = swpref.getpref('nmesh',[]);
nPatch0   = swpref.getpref('npatch',[]);

inpForm.fname  = {'range' 'mode'   'figure' 'color' 'linestyle' 'linewidth'};
inpForm.defval = {[]      'single' []       'auto'  '--'         1         };
inpForm.size   = {[-1 -2] [1 -3]   [1 1]    [1 -4]  [1 -5]       [1 1]     };
inpForm.soft   = {true    false    true     false   false        false     };

inpForm.fname  = [inpForm.fname  {'translate' 'zoom' 'tooltip' 'replace' 'legend'}];
inpForm.defval = [inpForm.defval { false       false true      true      true    }];
inpForm.size   = [inpForm.size   {[1 1]        [1 1] [1 1]     [1 1]     [1 1]   }];
inpForm.soft   = [inpForm.soft   {false        false false     false     false   }];

inpForm.fname  = [inpForm.fname  {'npatch' 'nmesh' 'unit' 'shift' 'copy'}];
inpForm.defval = [inpForm.defval {nPatch0  nMesh0  'lu'   [0;0;0] false }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]   [1 -6] [3 1]   [1 1] }];
inpForm.soft   = [inpForm.soft   {false    false   false  false   false }];

param = sw_readparam(inpForm, varargin{:});

if strcmp(param.color,'auto')
    param.color = [0 0 0];
end

if isempty(param.figure)
    hFigure  = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

% select range
if numel(param.range) == 6
    range = param.range;
elseif isempty(param.range)
    % get range from figure
    fRange = getappdata(hFigure,'range');
    if isempty(fRange)
        % fallback to default range
        range = [0 1;0 1;0 1];
    else
        % get plotting range and unit
        range       = fRange.range;
        param.unit  = fRange.unit;
    end
elseif numel(param.range) == 3
    % change range, if the number of unit cells are given
    range = [zeros(3,1) param.range(:)];
else
    error('plotcell:WrongInput','The given plotting range is invalid!');
end

switch param.mode
    case 'single'
        range = [0 1;0 1;0 1];
    case 'inside'
        range = [ceil(range(:,1)) floor(range(:,2))];
    case 'outside'
        range = [floor(range(:,1)) ceil(range(:,2))];
    otherwise
        error('plotcell:WrongInput','The given mode is invalid!');
end

% generate the unit cells
% TODO unit
Rz = cell(1,3);
[Rz{:}] = ndgrid(range(1,1):range(1,2),range(2,1):range(2,2),range(3,:));
Rz = reshape(permute(cat(4,Rz{:}),[4 1 2 3]),3,[],2);
Rx = cell(1,3);
[Rx{:}] = ndgrid(range(1,:),range(2,1):range(2,2),range(3,1):range(3,2));
Rx = reshape(permute(cat(4,Rx{:}),[4 2 3 1]),3,[],2);
Ry = cell(1,3);
[Ry{:}] = ndgrid(range(1,1):range(1,2),range(2,:),range(3,1):range(3,2));
Ry = reshape(permute(cat(4,Ry{:}),[4 1 3 2]),3,[],2);

pos = [Rx Ry Rz];

% basis vectors
BV = swplot.base(hFigure);

% shift positions
pos = bsxfun(@plus,pos,BV\param.shift);

% plot the cells already in base units, so no conversion needed
swplot.plot('type','line','position',pos,'figure',hFigure,...
    'linestyle',param.linestyle,'color',param.color,'name','cell',...
    'legend',false,'tooltip',false,'translate',param.translate,...
    'zoom',param.zoom,'replace',param.replace);

% save range
setappdata(hFigure,'range',struct('range',range,'unit',param.unit));

if nargout>0
    varargout{1} = hFigure;
end

if param.tooltip
    swplot.tooltip('on',hFigure);
end

end