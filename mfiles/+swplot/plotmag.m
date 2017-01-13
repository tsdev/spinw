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
% mode      String, defines the types of atoms to plot:
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
%           factor.
% R         Radius value of arrow body, default is 0.06.
% alpha     Head angle of the arrow in degree units, default is 30 degree.
% lHead     Length of the arrow head, default value is 0.5.
% aplane    Transparency (alpha value) of the circle, representing the
%           rotation plane of the moments, default is 0.07.
% centered  If true, the moment vector is centered on the atom, if false
%           the beggining of the spin vector is on the atom. Default is
%           true.
% nPatch    Number of points on the curve for the arrows, default
%           value is stored in swpref.getpref('npatch').
% tooltip   If true, the tooltips will be shown when clicking on atoms.
%           Default is true.
%
% Output:
%
% hFigure           Handle of the swplot figure.
%

% default values
fontSize0 = swpref.getpref('fontsize',[]);
nMesh0    = swpref.getpref('nmesh',[]);
nPatch0   = swpref.getpref('npatch',[]);
range0    = [0 1;0 1;0 1];

inpForm.fname  = {'range' 'legend' 'label' 'dtext' 'fontsize' 'radius0'};
inpForm.defval = {range0  true     true    0.1     fontSize0  0.3      };
inpForm.size   = {[-1 -2] [1 1]    [1 1]   [1 1]   [1 1]      [1 1]    };
inpForm.soft   = {false   false    false   false   false      false    };

inpForm.fname  = [inpForm.fname  {'radius' 'mode' 'color' 'nmesh' 'npatch'}];
inpForm.defval = [inpForm.defval {'auto'   'all'  'auto'  nMesh0  nPatch0 }];
inpForm.size   = [inpForm.size   {[1 -3]   [1 -4] [1 -5]  [1 1]   [1 1]   }];
inpForm.soft   = [inpForm.soft   {false    false  false   false   false   }];

inpForm.fname  = [inpForm.fname  {'hg'  'figure' 'obj' 'rangeunit' 'tooltip'}];
inpForm.defval = [inpForm.defval {true  []       []    'lu'        true     }];
inpForm.size   = [inpForm.size   {[1 1] [1 1]    [1 1] [1 -6]      [1 1]    }];
inpForm.soft   = [inpForm.soft   {false    true  true  false       false    }];

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

%lattice = obj.lattice;
% the basis vectors in columns.
BV = obj.basisvector;

% set figure title
set(hFigure,'Name', 'SpinW: Magnetic structure');

% change range, if the number of unit cells are given
if numel(param.range) == 3
    param.range = [ zeros(3,1) param.range(:)];
elseif numel(param.range) ~=6
    error('sw:plotatom:WrongInput','The given plotting range is invalid!');
end

range = param.range;

switch param.rangeunit
    case 'lu'
        rangelu = [floor(range(:,1)) ceil(range(:,2))];
    case 'xyz'
        % calculate the height of the paralellepiped of a unit cell
        hMax = 1./sqrt(sum(inv(BV).^2,2));
        
        % calculate the [111] position and find the cube that fits into the
        % parallelepiped
        hMax = min([sum(BV,2) hMax],[],2);
        
        % create range that includes the requested xyz box
        rangelu = bsxfun(@rdivide,range,hMax);
        rangelu = [floor(rangelu(:,1)) ceil(rangelu(:,2))];
        
    otherwise
        error('sw:plotatom:WrongInput','The given rangeunit string is invalid!');
end

% atom data
atom  = obj.atom;

% select atom type to plot
switch param.mode
    case 'all'
        % do nothing
    case 'mag'
        atom.r   = atom.r(:,atom.mag);
        atom.idx = atom.idx(1,atom.mag);
        atom.mag = atom.mag(1,atom.mag);
    case 'nonmag'
        atom.r   = atom.r(:,~atom.mag);
        atom.idx = atom.idx(1,~atom.mag);
        atom.mag = atom.mag(1,~atom.mag);
    otherwise
        error('sw:plotatom:WrongInput','The given mode string is invalid!');
end

nAtom = size(atom.r,2);

% generate positions from rangelu the inclusive range
pos = cell(1,3);
[pos{:}] = ndgrid(rangelu(1,1):rangelu(1,2),rangelu(2,1):rangelu(2,2),rangelu(3,1):rangelu(3,2));
pos = bsxfun(@plus,reshape(cat(4,pos{:}),[],3)',permute(atom.r,[1 3 2]));

nCell = size(pos,2);

% keep track of types of atoms
aIdx = repmat(atom.idx,[nCell 1]);

pos  = reshape(pos,3,[]);
aIdx = reshape(aIdx,3,[]);

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
    warning('plotatom:EmptyPlot','There are no atoms in the plotting range!')
    return
end

pos  = pos(:,pIdx);
aIdx = aIdx(pIdx);

% color
if strcmp(param.color,'auto')
    color = double(obj.unit_cell.color(:,aIdx));
else
    color = swplot.color(param.color);
end

% radius
if strcmp(param.radius,'auto')
    radius = sw_atomdata(obj.unit_cell.label,'radius');
    radius = radius(aIdx)*param.radius0;
else
    radius = param.radius0;
end

% prepare labels
atom.name = obj.unit_cell.label(atom.idx);
% plot only the first word of every label
for ii = 1:nAtom
    labelTemp = strword(atom.name{ii},[1 2],true);
    label1 = labelTemp{1};
    %label2 = labelTemp{2};
    atom.text{ii}  = [label1 '(' num2str(atom.idx(ii)) ')_' num2str(ii)];
    %atom.ttip{ii}  = [label2 ' atom (' label1 ')' newline 'Unit cell:' newline];
end

% plot label on atoms
if param.label
    dtext = inv(BV)*repmat(radius+param.dtext,[3 1]); %#ok<MINV>
    text = repmat(atom.text,[nCell 1]);
    text = text(pIdx);
    
    swplot.plot('type','text','name','atom_label','position',pos+dtext,...
        'text',text,'figure',hFigure,'legend',false,'fontsize',param.fontsize,'tooltip',false);
end

% legend label
lLabel = repmat(atom.name,[nCell 1]);
lLabel = lLabel(pIdx);

if param.legend
    
    % create specific legend per atom type
    lDat = getappdata(hFigure,'legend');
    % append color
    lDat.color = [lDat.color double(obj.unit_cell.color)/255];
    % append type
    lDat.type = [lDat.type 3*ones(1,obj.natom)];
    % append text
    if ~isempty(lDat.text)
        lDat.text = [lDat.text(:)' obj.unit_cell.label];
    else
        lDat.text = obj.unit_cell.label;
    end
    
    setappdata(hFigure,'legend',lDat);
    swplot.legend('on',hFigure);
end

% plot the atoms, text generated automatically
swplot.plot('type','ellipsoid','name','atom','position',pos,'R',radius,...
    'figure',hFigure,'color',color,'text','','legend',false,'label',lLabel,...
    'nmesh',param.nmesh,'tooltip',false);

if nargout > 0
    varargout{1} = hFigure;
end

if param.tooltip
    swplot.tooltip('on',hFigure);
end

end