function varargout = plotchem(varargin)
% plots polyhedra or chemical bonds
%
% SWPLOT.PLOTCHEM('option1', value1, ...)
%
% hFigure = SWPLOT.PLOTCHEM(...)
%
% The function polyhedra around selected  atoms, or chemical bonds between
% atoms an swplot figure.
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
% unit      Unit in which the range is defined. It can be the following
%           string:
%               'lu'        Lattice units (default).
%               'xyz'       Cartesian coordinate system in Angstrom units.
% mode      Selects the type of the plot:
%               'poly'      Draws polyhedra around the center atoms
%                           (default).
%               'bond'      Draws bonds between given atoms.
% atom1     Indices of atoms stored in spinw.unit_cell for the center atom
%           of the polyhedra or the first atom of the bonds. Can be also a
%           string that identifies the atoms by their labels.
% atom2     Indices or label of the atoms stored in spinw.unit_cell. It
%           determines the vertices of the polyhedra or gives the second
%           atom of a bond.
% limit     Can be a single number which will restrict the number of
%           nearest negihbours of atom1 to connect. Can be also a vector
%           that defines bonds/polyhedra between atoms that are within the
%           given distance range stored as a row vector [dmin dmax].
%           Default is 6 to plot octahedra around atom1.
% alpha     Transparency of the plotted surfaces between 0 and 1 (1 for
%           opaque, 0 for transparent). Default value is 1 for bonds and
%           0.5 for polyhedron.
% color     Surface color of the objects. Default is 'auto', when they are
%           set to the color of atom1. [R G B] will fix the color of all
%           bonds to a uniform one, can be also arbitrary color name (see
%           swplot.color() function). Can be also 'none', when no faces
%           will be shown.
% color2    Color of the edges of the polyhedra (unused for bonds), default
%           value is 'auto' when the edge gets the same color as the faces.
%           'none' will remove the edges.
% radius0   Radius of the cylinder, default value is 0.03.
% figure    Handle of the swplot figure. Default is the selected figure.
% legend    Whether to add the plot to the legend, default is true.
% nPatch    Number of points on the curve for the cylinder, default
%           value is stored in swpref.getpref('npatch').
% tooltip   If true, the tooltips will be shown when clicking on atoms.
%           Default is true.
% shift     Column vector with 3 elements, all atomic positions will be
%           shifted by the given value. Default value is [0;0;0].
% replace   Replace previous atom plot if true. Default is true.
% translate If true, all plot objects will be translated to the figure
%           center. Default is true.
% zoom      If true, figure will be automatically zoomed to the ideal size.
%           Default is true.
%
% Output:
%
% hFigure           Handle of the swplot figure.
%
% The name of the objects that are created called 'chem'. To find the
% handles and the stored data on these objects, use e.g.
%
%   sObject = swplot.findobj(hFigure,'name','chem')
%


% default values
nMesh0    = swpref.getpref('nmesh',[]);
nPatch0   = swpref.getpref('npatch',[]);
range0    = [0 1;0 1;0 1];

inpForm.fname  = {'range' 'legend' 'label' 'unit' 'mode' 'atom1' 'atom2'};
inpForm.defval = {range0  true     true    'lu'   'poly' 1       2      };
inpForm.size   = {[-1 -2] [1 1]    [1 1]   [1 -3] [1 -4] [1 -5]  [1 -6] };
inpForm.soft   = {false   false    false   false  false  false   false  };

inpForm.fname  = [inpForm.fname  {'limit' 'alpha' 'color' 'nmesh' 'npatch'}];
inpForm.defval = [inpForm.defval {6       1       'auto'  nMesh0  nPatch0 }];
inpForm.size   = [inpForm.size   {[1 -7]  [1 1]   [1 -8]  [1 1]   [1 1]   }];
inpForm.soft   = [inpForm.soft   {false   false   false   false   false   }];

inpForm.fname  = [inpForm.fname  {'figure' 'obj' 'color2' 'tooltip' 'radius0'}];
inpForm.defval = [inpForm.defval {[]       []    'auto'   true      0.03     }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1] [1 -9]  [1 1]      [1 1]    }];
inpForm.soft   = [inpForm.soft   {true     true  false   false      false    }];

inpForm.fname  = [inpForm.fname  {'shift' 'replace' 'translate' 'zoom'}];
inpForm.defval = [inpForm.defval {[0;0;0] true      true         true }];
inpForm.size   = [inpForm.size   {[3 1]   [1 1]     [1 1]        [1 1]}];
inpForm.soft   = [inpForm.soft   {false   false     false        false}];

param = sw_readparam(inpForm, varargin{:});

if isempty(param.figure)
    hFigure  = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

if isempty(param.obj) && ~isappdata(hFigure,'obj')
    warning('plotchem:WrongInput','No SpinW object to plot!');
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

% the basis vectors in columns.
BV = obj.basisvector;

% set figure title
set(hFigure,'Name', 'SpinW: Chemical bonds');

% change range, if the number of unit cells are given
if numel(param.range) == 3
    param.range = [ zeros(3,1) param.range(:)];
elseif numel(param.range) ~=6
    error('plotchem:WrongInput','The given plotting range is invalid!');
end

range = param.range;

switch param.unit
    case 'lu'
        rangelu = [floor(range(:,1)) ceil(range(:,2))];
    case 'xyz'
        % corners of the box
        corners = BV\[range(1,[1 2 2 1 1 2 2 1]);range(2,[1 1 2 2 1 1 2 2]);range(3,[1 1 1 1 2 2 2 2])];
        rangelu = [min(corners,[],2) max(corners,[],2)];
        rangelu = [floor(rangelu(:,1)) ceil(rangelu(:,2))];
        
    otherwise
        error('plotchem:WrongInput','The given unit string is invalid!');
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
        error('plotchem:WrongInput','The given mode string is invalid!');
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
%aIdx = reshape(aIdx,3,[]);

% cut out the atoms that are out of range
switch param.unit
    case 'lu'
        % L>= lower range, L<= upper range
        pIdx = all(bsxfun(@ge,pos,range(:,1)) & bsxfun(@le,pos,range(:,2)),1);
    case 'xyz'
        % convert to xyz
        posxyz = BV*pos;
        pIdx = all(bsxfun(@ge,posxyz,range(:,1)) & bsxfun(@le,posxyz,range(:,2)),1);
end

if ~any(pIdx)
    warning('plotchem:EmptyPlot','There are no atoms in the plotting range!')
    return
end

pos  = pos(:,pIdx);
aIdx = aIdx(pIdx);
aIdx = aIdx(:)';

% color
if strcmp(param.color,'auto')
    color = double(obj.unit_cell.color(:,aIdx));
else
    color = swplot.color(param.color);
end

% radius
switch param.radius
    case 'auto'
        radius = sw_atomdata(obj.unit_cell.label,'radius');
        radius = radius(aIdx)*param.radius0;
    case 'fix'
        radius = param.radius0;
    otherwise
        error('plotchem:WrongInput','The given radius option is invalid!');
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
pos = bsxfun(@plus,pos,param.shift);

% plot label on atoms
if param.label
    dtext = inv(BV)*repmat(radius+param.dtext,[3 1]); %#ok<MINV>
    text = repmat(atom.text,[nCell 1]);
    text = text(pIdx);
    
    swplot.plot('type','text','name','atom_label','position',pos+dtext,...
        'text',text,'figure',hFigure,'legend',false,'translate',false,'zoom',false,...
        'fontsize',param.fontsize,'tooltip',false,'replace',param.replace);
end

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
    setappdata(hFigure,'legend',lDat);
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
    'nmesh',param.nmesh,'tooltip',false,'data',posDat,'replace',param.replace,...
    'translate',param.translate,'zoom',param.zoom);

if nargout > 0
    varargout{1} = hFigure;
end

% save range
setappdata(hFigure,'range',struct('range',param.range,'unit',param.unit));

if param.tooltip
    swplot.tooltip('on',hFigure);
end

end