function varargout = plotcell(varargin)
% plots the edges of unit cells on swplot figure
%
% SWPLOT.PLOTCELL('Option1',Value1,...)
%
% hFigure = SWPLOT.PLOTCELL('Option1',Value1,...)
%
% Options:
%
% range     Plotting range of the lattice parameters in lattice units,
%           dimensions are [3 2]. For example to plot the first unit cell,
%           use: [0 1;0 1;0 1]. Also the number unit cells can be given
%           along the a, b and c directions: [2 1 2], that is equivalent to
%           [0 2;0 1;0 2]. Default is the single unit cell.
% mode      String, determined how the cells are plotted:
%               'single'    A single unit cell is plotted at the origin.
%               'inside'    Unit cells are plotted inside the given
%                           range.
%               'outside'   Unit cells are plotted inclusive the given
%                               range.
% figure    Handle of the swplot figure. Default is the selected figure.
% color     Color of the lines as a string or row vector with 3 elements, 
%           default value is black ('auto').
% lineStyle Line style of the cell, default is '--'.
% lineWdith Line width of the cell, default is 1.
% translate If true, all plot objects will be translated to the figure
%           center. Default is true.
% zoom      If true, figure will be automatically zoomed to the ideal size.
%           Default is true.
%

% default values
range0    = [0 1;0 1;0 1];
nMesh0    = swpref.getpref('nmesh',[]);
nPatch0   = swpref.getpref('npatch',[]);

inpForm.fname  = {'range' 'mode'   'figure' 'color' 'linestyle' 'linewidth'};
inpForm.defval = {range0  'single' []       'auto'  '--'         1         };
inpForm.size   = {[-1 -2] [1 -3]   [1 1]    [1 -4]  [1 -5]       [1 1]     };
inpForm.soft   = {false   false    true     false   false        false     };

inpForm.fname  = [inpForm.fname  {'translate' 'zoom' 'tooltip' 'replace'}];
inpForm.defval = [inpForm.defval { true         true true      true     }];
inpForm.size   = [inpForm.size   {[1 1]        [1 1] [1 1]     [1 1]    }];
inpForm.soft   = [inpForm.soft   {false        false false     false    }];

inpForm.fname  = [inpForm.fname  {'npatch' 'nmesh' 'unit'}];
inpForm.defval = [inpForm.defval {nPatch0  nMesh0  'lu'  }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]   [1 -6]}];
inpForm.soft   = [inpForm.soft   {false    false   false }];

param = sw_readparam(inpForm, varargin{:});

if strcmp(param.color,'auto')
    param.color = [0 0 0];
end

if isempty(param.figure)
    hFigure  = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

% change range, if the number of unit cells are given
if numel(param.range) == 3
    param.range = [ zeros(3,1) param.range(:)];
elseif numel(param.range) ~=6
    error('sw:cell:WrongInput','The given plotting range is invalid!');
end

switch param.mode
    case 'single'
        range = range0;
    case 'inside'
        range = [ceil(param.range(:,1)) floor(param.range(:,2))];
    case 'outside'
        range = [floor(param.range(:,1)) ceil(param.range(:,2))];
    otherwise
        error('sw:cell:WrongInput','The given mode is invalid!');
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

R = [Rx Ry Rz];

% plot the cells already in base units, so no conversion needed
swplot.plot('type','line','position',R,'figure',hFigure,...
    'linestyle',param.linestyle,'color',param.color,'name','cell',...
    'legend',false,'tooltip',false,'translate',param.translate,...
    'zoom',param.zoom,'replace',param.replace);

% save range
setappdata(hFigure,'range',struct('range',param.range,'unit',param.unit));

if nargout>0
    varargout{1} = hFigure;
end

if param.tooltip
    swplot.tooltip('on',hFigure);
end

end