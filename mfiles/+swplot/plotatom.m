function varargout = plotatom(varargin)
% plots crystal structure
%
% SWPLOT.PLOTATOM('option1', value1, ...)
%
% hFigure = SWPLOT.PLOTATOM(...)
%
% The function plots the crystal structure of a SpinW object onto an swplot
% figure.
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
%               'all'       Plot all atoms, default value.
%               'mag'       Plot magnetic atoms only.
%               'nonmag'    Plot non-magnetic atoms only.
% figure    Handle of the swplot figure. Default is the selected figure.
% legend    Whether to add the plot to the legend, default is true.
% label     Whether to plot labels for atoms, default is true.
% dText     Distance between item and its text label, default is 0.1
%           Angstrom.
% fontSize  Font size of the atom labels in pt, default is stored in
%           swpref.getpref('fontsize').
% radius0   Constant atom radius, default value is 0.3 Angstrom.
% radius    Defines the atom radius:
%               numerical value     Sets the radius of all atoms to the
%                                   given value. Default is radius0.
%               'auto'              use radius data from database based on
%                                   the atom label multiplied by radius0
%                                   value.
% color     Color of the atoms:
%               'auto'      All atom gets the color stored in obj.unit_cell.
%               'colorname' All atoms will have the same color.
%               [R G B]     RGB code of the color that fix the color of all
%                           atoms.
% nMesh     Resolution of the ellipse surface mesh. Integer number that is
%           used to generate an icosahedron mesh with #mesh number of
%           additional triangulation, default value is stored in
%           swpref.getpref('nmesh')
% tooltip   If true, the tooltips will be shown when clicking on atoms.
%           Default is true.
% shift     Column vector with 3 elements, all atomic positions will be
%           shifted by the given value. Default value is [0;0;0].
% replace   Replace previous atom plot if true. Default is true.
% hg        Whether to use hgtransform (nice rotation with the mouse) or
%           default Matlab rotation of 3D objects. Default is true.
%
% Output:
%
% For plotting:
% hFigure           Handle of the swplot figure.
%
% The labels on the atom has the following meaning:
%
% AtomLabel(idx1)_idx2:
%   idx1    The index of the atom in obj.unit_cell
%   idx2    The index of the atom in obj.atom
%
% To access the graphical objects, use the function:
%   handleList = sw_getobject(tagName,fHandle);
% This command returns all objects with tagName string in their 'Tag'
% field, in figure with fHandle.
%
% List of Tag values:
% unitCell         Unit cell.
% arrow            Arrows of the lattice vectors.
% arrowText        Text of the lattice vectors.
% atom_            Spheres of the atoms, for example 'atom_Cr' for atoms
%                  with the label 'Cr'.
% aniso_           Anisotropy ellipsoid with atom names, for example
%                  'aniso_Cr', symbolising the anisotropy matrix.
% atomText_        Label of atoms, for example 'atomText_Cr'.
% coupling_        Cylinder of couplings, for example 'coupling_J1' for the
%                  couplings with label 'J1'.
% spinPlane        Circle showing the plane of the spins.
% legend           Legend.
%
% Example:
%
% sq = sw_model('squareAF',1);
% plot(sq);
%
% The example plots the structure, couplings and magnetic structure of the
% square lattice antiferromagnet.
%
% See also SW_DRAW, SW_GETOBJECT, SW_ADDOBJECT.
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

inpForm.fname  = [inpForm.fname  {'shift' 'replace' }];
inpForm.defval = [inpForm.defval {[0;0;0] true      }];
inpForm.size   = [inpForm.size   {[3 1]   [1 1]     }];
inpForm.soft   = [inpForm.soft   {false   false     }];

param = sw_readparam(inpForm, varargin{:});

if isempty(param.figure)
    hFigure  = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

if isempty(param.obj) && ~isappdata(hFigure,'obj')
    warning('plotatom:WrongInput','No SpinW object to plot!');
    return
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
set(hFigure,'Name', 'SpinW: Crystal structure');

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
        hMax = 1./sqrt(sum(inv(BV').^2,2));
        
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

% save atom coordinates into data
posDat = mat2cell(pos,3,ones(1,size(pos,2)));

% shift positions
pos = pos+param.shift;

% plot label on atoms
if param.label
    dtext = inv(BV)*repmat(radius+param.dtext,[3 1]); %#ok<MINV>
    text = repmat(atom.text,[nCell 1]);
    text = text(pIdx);
    
    swplot.plot('type','text','name','atom_label','position',pos+dtext,...
        'text',text,'figure',hFigure,'legend',false,...
        'fontsize',param.fontsize,'tooltip',false,'replace',param.replace);
end

% % tooltip
% ttip = repmat(atom.ttip,[nCell 1]);
% ttip = ttip(pIdx);
% % add cell index and position
% for ii = 1:numel(ttip)
%     posi = pos(:,ii)';
%     ttip{ii} = [ttip{ii} sprintf('[%d, %d, %d]',floor(posi)) newline 'Atomic position' newline sprintf('[%5.3f, %5.3f, %5.3f]',posi-floor(posi))];
% end

% legend label
lLabel = repmat(atom.name,[nCell 1]);
lLabel = lLabel(pIdx);

% legend data
lDat = getappdata(hFigure,'legend');

if param.replace
    % remove old legend entries
    lIdx = ~ismember(lDat.name,'atom');
    lDat.color = lDat.color(:,lIdx);
    lDat.type  = lDat.type(:,lIdx);
    lDat.name  = lDat.name(:,lIdx);
    lDat.text  = lDat.text(:,lIdx);
    % redraw legend
    swplot.legend('refresh',hFigure);
end

if param.legend
    % append color
    lDat.color = [lDat.color double(obj.unit_cell.color)/255];
    % append type
    lDat.type = [lDat.type 3*ones(1,obj.natom)];
    % append name
    lDat.name = [lDat.name repmat({'atom'},1,obj.natom)];
    % append text
    lDat.text = [lDat.text obj.unit_cell.label];
    
    setappdata(hFigure,'legend',lDat);
    swplot.legend('on',hFigure);
end

% plot the atoms, text generated automatically
swplot.plot('type','ellipsoid','name','atom','position',pos,'R',radius,...
    'figure',hFigure,'color',color,'text','','legend',false,'label',lLabel,...
    'nmesh',param.nmesh,'tooltip',false,'data',posDat,'replace',param.replace);

if nargout > 0
    varargout{1} = hFigure;
end

if param.tooltip
    swplot.tooltip('on',hFigure);
end

end