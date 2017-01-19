function varargout = plotmag(varargin)
% plots magnetic structure
%
% SWPLOT.PLOTMAG('option1', value1, ...)
%
% hFigure = SWPLOT.PLOTMAG(...)
%
% The function plots the magnetic structure of a SpinW object onto an
% swplot figure.
%
% Input:
%
% Options:
%
% obj       SpinW object.
% range     Plotting range of the lattice parameters in lattice units,
%           dimensions are [3 2]. For example to plot the first unit cell,
%           use: [0 1;0 1;0 1]. Also the number unit cells can be given
%           along the a, b and c directions: [2 1 2], that is equivalent to
%           [0 2;0 1;0 2]. Default is the single unit cell.
% rangeunit Unit in which the range is defined. It can be the following
%           string:
%               'lu'        Lattice units (default).
%               'xyz'       Cartesian coordinate system in Angstrom units.
% mode      String, defines the way the magnetic moments are plotted:
%               'all'       Plot both the rotation plane of incommensurate
%                           magnetic structures and the moment directions.
%               'circle'    Plot only the rotation plane of incommensurate
%                           magnetic structures.
%               'arrow'     Plots only the moment directions.
% figure    Handle of the swplot figure. Default is the selected figure.
% legend    Whether to add the plot to the legend, default is true.
% label     Whether to plot labels for atoms, default is true.
% dText     Distance between item and its text label, default is 0.1
%           Angstrom.
% fontSize  Font size of the atom labels in pt, default value is stored in
%           swpref.getpref('fontsize').
% color     Color of the magnetic moments:
%               'auto'      All moments get the same color as the magnetic
%                           atom.
%               'colorname' All moments will have the same color.
%               [R G B]     RGB code of the color.
% scale     Scaling factor for the lenght of the magnetic moments.
% normalize If true, all moment length will be normalized to the scale
%           factor, default is false.
% radius0   Radius value of arrow body, default is 0.06.
% ang       Angle of the arrow head in degree units, default is 30 degree.
% lHead     Length of the arrow head, default value is 0.5.
% alpha     Transparency (alpha value) of the circle, representing the
%           rotation plane of the moments, default is 0.07.
% centered  If true, the moment vector is centered on the atom, if false
%           the beggining of the spin vector is on the atom. Default is
%           true.
% nPatch    Number of points on the curve for the arrows, default
%           value is stored in swpref.getpref('npatch').
% tooltip   If true, the tooltips will be shown when clicking on atoms.
%           Default is true.
% shift     Column vector with 3 elements, all vectors will be
%           shifted by the given value. Default value is [0;0;0].
% replace   Replace previous magnetic moment plot if true. Default is true.
% translate If true, all plot objects will be translated to the figure
%           center. Default is true.
% zoom      If true, figure will be automatically zoomed to the ideal size.
%           Default is true.
%
% Output:
%
% hFigure           Handle of the swplot figure.
%
% The name of the objects that are created called 'mag'.
% To find the handles and the stored data on these objects, use e.g.
%
%   sObject = swplot.findobj(hFigure,'name','mag')
%


% default values
fontSize0 = swpref.getpref('fontsize',[]);
nMesh0    = swpref.getpref('nmesh',[]);
nPatch0   = swpref.getpref('npatch',[]);
range0    = [0 1;0 1;0 1];

inpForm.fname  = {'range' 'legend' 'label' 'dtext' 'fontsize' 'radius0' };
inpForm.defval = {range0  true     true    0.1     fontSize0  0.06      };
inpForm.size   = {[-1 -2] [1 1]    [1 1]   [1 1]   [1 1]      [1 1]     };
inpForm.soft   = {false   false    false   false   false      false     };

inpForm.fname  = [inpForm.fname  {'mode' 'color' 'nmesh' 'npatch' 'ang' }];
inpForm.defval = [inpForm.defval {'all'  'auto'  nMesh0  nPatch0  30    }];
inpForm.size   = [inpForm.size   {[1 -4] [1 -5]  [1 1]   [1 1]    [1 1] }];
inpForm.soft   = [inpForm.soft   {false  false   false   false    false }];

inpForm.fname  = [inpForm.fname  {'figure' 'obj' 'rangeunit' 'tooltip'}];
inpForm.defval = [inpForm.defval {[]       []    'lu'        true     }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1] [1 -6]      [1 1]    }];
inpForm.soft   = [inpForm.soft   {true     true  false       false    }];

inpForm.fname  = [inpForm.fname  {'shift' 'replace' 'scale' 'normalize' }];
inpForm.defval = [inpForm.defval {[0;0;0] true      1        false      }];
inpForm.size   = [inpForm.size   {[3 1]   [1 1]     [1 1]    [1 1]      }];
inpForm.soft   = [inpForm.soft   {false   false     false    false      }];

inpForm.fname  = [inpForm.fname  {'lHead' 'alpha' 'centered' 'translate' 'zoom'}];
inpForm.defval = [inpForm.defval {0.5      0.07   false      true         true }];
inpForm.size   = [inpForm.size   {[1 1]   [1 1]   [1 1]      [1 1]        [1 1]}];
inpForm.soft   = [inpForm.soft   {false   false   false      false        false}];

% TODO param.scale

param = sw_readparam(inpForm, varargin{:});

if isempty(param.figure)
    hFigure  = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

% takes care of spinw object saved/loaded in/from figure
if isempty(param.obj)
    obj = getappdata(hFigure,'obj');
else
    setappdata(hFigure,'obj',copy(param.obj));
    obj = param.obj;
    setappdata(hFigure,'base',obj.basisvector);
end

% the basis vectors in columns.
BV = obj.basisvector;

% set figure title
set(hFigure,'Name', 'SpinW: Magnetic structure');

% change range, if the number of unit cells are given
if numel(param.range) == 3
    param.range = [ zeros(3,1) param.range(:)];
elseif numel(param.range) ~=6
    error('plotmag:WrongInput','The given plotting range is invalid!');
end

range = param.range;

switch param.rangeunit
    case 'lu'
        rangelu = [floor(range(:,1)) ceil(range(:,2))];
    case 'xyz'
        % corners of the box
        corners = BV\[range(1,[1 2 2 1 1 2 2 1]);range(2,[1 1 2 2 1 1 2 2]);range(3,[1 1 1 1 2 2 2 2])];
        rangelu = [min(corners,[],2) max(corners,[],2)];
        rangelu = [floor(rangelu(:,1)) ceil(rangelu(:,2))];
    otherwise
        error('plotmag:WrongInput','The given rangeunit string is invalid!');
end

% magnetic supercell for plotting
nExtPlot = diff(rangelu,1,2)'+1;
luOrigin = rangelu(:,1)';
% magnetic moment data
magStr   = obj.magstr('nExt',nExtPlot,'origin',luOrigin);
M        = magStr.S;

% magnetic atom position data
mAtom    = obj.matom;
%nMagAtom = size(mAtom.r,2);

% atom positions in the plotting range
mAtomExt = sw_extendlattice(nExtPlot,mAtom);

pos = bsxfun(@plus,bsxfun(@times,mAtomExt.RRext,nExtPlot'),luOrigin');

% select atom type to plot
switch param.mode
    case 'all'
        % do nothing
    case 'circle'
    case 'arrow'
    otherwise
        error('plotmag:WrongInput','The given mode string is invalid!');
end

% number of unit cells
nCell = prod(nExtPlot);

% keep track of types of atoms
aIdx = repmat(mAtom.idx,[nCell 1])';
pos  = reshape(pos,3,[]);

% cut out the atoms that are out of range
switch param.rangeunit
    case 'lu'
        % L>= lower range, L<= upper range
        pIdx = all(bsxfun(@ge,pos,range(:,1)) & bsxfun(@le,pos,range(:,2)),1);
    case 'xyz'
        % convert to xyz
        posxyz = BV*pos;
        pIdx = all(bsxfun(@ge,posxyz,range(:,1)) & bsxfun(@le,posxyz,range(:,2)),1);
end

if ~any(pIdx)
    warning('plotatom:EmptyPlot','There are no magnetic moments in the plotting range!')
    return
end

M    = M(:,pIdx);
pos  = pos(:,pIdx);
aIdx = aIdx(pIdx);

% normalization
if param.normalize
    % normalize moments
    M = bsxfun(@rdivide,M,sqrt(sum(M.^2,1)));
end

% save magnetic moment vector values into data before rescaling
MDat = mat2cell(M,3,ones(1,size(M,2)));

% scale moments
M = M*param.scale;

% convert to lu units
Mlu = BV\M;

if param.centered
    % center on atoms
    vpos = cat(3,pos-Mlu/2,pos+Mlu/2);
else
    vpos = cat(3,pos,pos+Mlu);
end

% color
if strcmp(param.color,'auto')
    color = double(obj.unit_cell.color(:,aIdx));
else
    color = swplot.color(param.color);
end

% shift positions
vpos = bsxfun(@plus,vpos,param.shift);

% prepare legend labels
mAtom.name = obj.unit_cell.label(mAtom.idx);
lLabel = repmat(mAtom.name,[nCell 1]);
lLabel = lLabel(pIdx);


% plot moment vectors
swplot.plot('type','arrow','name','mag','position',vpos,'text','',...
    'figure',hFigure,'legend',false,'color',color,'R',param.radius0,...
    'fontsize',param.fontsize,'tooltip',false,'replace',param.replace,...
    'data',MDat,'label',lLabel,'nmesh',param.nmesh,'ang',param.ang,...
    'lHead',param.lHead,'translate',param.translate,'zoom',param.zoom);

% save range
setappdata(hFigure,'range',param.range);

if nargout > 0
    varargout{1} = hFigure;
end

if param.tooltip
    swplot.tooltip('on',hFigure);
end

end