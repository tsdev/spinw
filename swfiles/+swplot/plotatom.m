function varargout = plotatom(varargin)
% plots crystal structure
% 
% ### Syntax
% 
% `swplot.plotatom(Name,Value)`
% 
% `hFigure = swplot.plotatom(Name,Value)`
%
% ### Description
% 
% `swplot.plotatom(Name,Value)` plots the crystal structure of a SpinW
% object onto an swplot figure where each atom is shown as a sphere. It can
% display text labels, control the radius and color of the spheres.
%  
% `hFigure = swplot.plotatom(Name,Value)` returns the handle of the swplot
% figure.
%
% ### Name-Value Pair Arguments
% 
% `'mode'`
% : String that defines the types of atoms to plot:
%   * `'all'`       plot all atoms (default),
%   * `'mag'`       plot magnetic atoms only,
%   * `'nonmag'`    plot non-magnetic atoms only.
% 
% `'label'`
% : Whether to plot the labels of the atoms, default value is `false`.
% 
% `'dText'`
% : Distance between the atom and its text label, default value is 0.1
%   \\ang.
% 
% `'fontSize'`
% : Font size of the atom labels in pt, default value is read using
%   `swpref.getpref('fontsize')`.
% 
% `'radius'`
% : Defines the atom radius, one of the following strings:
%   * `'fix'`       Sets the radius of all atoms to the value
%                   of the `radius0` parameter,
%   * `'auto'`      determine the atom radius based on the atom
%                   label (e.g. the radius of Cr atom is
%                   `sw_atomdata('Cr','radius')`) and multiply the value by
%                   the `radius0` parameter.
%
% `'radius0'`
% : Constant atom radius, default value is 0.3 \\ang.
% 
% `'color'`
% : Color of the atoms, one of the following values:
%   * `'auto'`      all atom gets the color stored in [spinw.unit_cell],
%   * `'colorname'` all atoms will have the same color defined by the
%                   string, e.g. `'red'`,
%   * `[R G B]`     all atoms will have the same color defined by the RGB
%                   code.
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
% The name of the objects that are created are `'atom'` and `'atom_label'`.
% To find the handles and the corresponding data on these objects, use e.g.
% sObject = swplot.findobj(hFigure,'name','atom')`.
%

% default values
fontSize0 = swpref.getpref('fontsize',[]);
nMesh0    = swpref.getpref('nmesh',[]);
nPatch0   = swpref.getpref('npatch',[]);

inpForm.fname  = {'range' 'legend' 'label' 'dtext' 'fontsize' 'radius0'};
inpForm.defval = {[]      true     false    0.1     fontSize0  0.3     };
inpForm.size   = {[-1 -2] [1 1]    [1 1]   [1 1]   [1 1]      [1 1]    };
inpForm.soft   = {true    false    false   false   false      false    };

inpForm.fname  = [inpForm.fname  {'radius' 'mode' 'color' 'nmesh' 'npatch'}];
inpForm.defval = [inpForm.defval {'auto'   'all'  'auto'  nMesh0  nPatch0 }];
inpForm.size   = [inpForm.size   {[1 -3]   [1 -4] [1 -5]  [1 1]   [1 1]   }];
inpForm.soft   = [inpForm.soft   {false    false  false   false   false   }];

inpForm.fname  = [inpForm.fname  {'figure' 'obj' 'unit'  'tooltip' 'copy'}];
inpForm.defval = [inpForm.defval {[]       []    'lu'    true      false }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1] [1 -6]  [1 1]     [1 1] }];
inpForm.soft   = [inpForm.soft   {true     true  false   false     false }];

inpForm.fname  = [inpForm.fname  {'shift' 'replace' 'translate' 'zoom'}];
inpForm.defval = [inpForm.defval {[0;0;0] true      false        false}];
inpForm.size   = [inpForm.size   {[3 1]   [1 1]     [1 1]        [1 1]}];
inpForm.soft   = [inpForm.soft   {false   false     false        false}];

param = sw_readparam(inpForm, varargin{:});

if isempty(param.figure)
    hFigure  = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

% take care of output
if nargout > 0
    varargout{1} = hFigure;
end

if isempty(param.obj) && ~isappdata(hFigure,'obj')
    warning('plotatom:WrongInput','No SpinW object to plot!');
    return
end

% takes care of spinw object saved/loaded in/from figure
if isempty(param.obj)
    obj = getappdata(hFigure,'obj');
else
    if param.copy
        setappdata(hFigure,'obj',copy(param.obj));
    else
        setappdata(hFigure,'obj',param.obj);
    end
    obj = param.obj;
    setappdata(hFigure,'base',obj.basisvector);
end

%lattice = obj.lattice;
% the basis vectors in columns.
BV = obj.basisvector;

% set figure title
set(hFigure,'Name', 'SpinW: Crystal structure');

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
    error('plotatom:WrongInput','The given plotting range is invalid!');
end

switch param.unit
    case 'lu'
        rangelu = [floor(range(:,1)) ceil(range(:,2))];
    case 'xyz'
        % corners of the box
        corners = BV\[range(1,[1 2 2 1 1 2 2 1]);range(2,[1 1 2 2 1 1 2 2]);range(3,[1 1 1 1 2 2 2 2])];
        rangelu = [min(corners,[],2) max(corners,[],2)];
        rangelu = [floor(rangelu(:,1)) ceil(rangelu(:,2))];
        
    otherwise
        error('plotatom:WrongInput','The given unit string is invalid!');
end

% atom data
atom  = obj.atom;

% select atom type to plot
switch param.mode
    case 'all'
        % do nothing
        lSelect = true(size(obj.unit_cell.S));
        a2Idx = 1:numel(atom.idx);
    case 'mag'
        a2Idx    = find(atom.mag);
        atom.r   = atom.r(:,atom.mag);
        atom.idx = atom.idx(1,atom.mag);
        atom.mag = atom.mag(1,atom.mag);
        lSelect  = obj.unit_cell.S>0;
    case 'nonmag'
        a2Idx    = find(~atom.mag);
        atom.r   = atom.r(:,~atom.mag);
        atom.idx = atom.idx(1,~atom.mag);
        atom.mag = atom.mag(1,~atom.mag);
        lSelect  = obj.unit_cell.S==0;
    otherwise
        error('plotatom:WrongInput','The given mode string is invalid!');
end

nAtom = size(atom.r,2);

% generate positions from rangelu the inclusive range
pos = cell(1,3);
[pos{:}] = ndgrid(rangelu(1,1):rangelu(1,2),rangelu(2,1):rangelu(2,2),rangelu(3,1):rangelu(3,2));
pos = bsxfun(@plus,reshape(cat(4,pos{:}),[],3)',permute(atom.r,[1 3 2]));

nCell = size(pos,2);

% keep track of types of atoms (index in spinw.unit_cell)
aIdx = repmat(atom.idx,[nCell 1]);
% keep track of atom index in spinw.atom list
a2Idx = repmat(a2Idx,[nCell 1]);

pos  = reshape(pos,3,[]);
%aIdx = reshape(aIdx,3,[]);

% cut out the atoms that are out of range
switch param.unit
    case 'lu'
        % L>= lower range, L<= upper range
        pIdx = all(bsxfun(@ge,pos,range(:,1)-10*eps) & bsxfun(@le,pos,range(:,2)+10*eps),1);
    case 'xyz'
        % convert to xyz
        posxyz = BV*pos;
        pIdx = all(bsxfun(@ge,posxyz,range(:,1)-10*eps) & bsxfun(@le,posxyz,range(:,2)+10*eps),1);
end

if ~any(pIdx)
    warning('plotatom:EmptyPlot','There are no atoms in the plotting range!')
    return
end

pos   = pos(:,pIdx);
aIdx  = aIdx(pIdx);
aIdx  = aIdx(:)';
a2Idx = a2Idx(pIdx);
a2Idx = a2Idx(:)';

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
        error('plotatom:WrongInput','The given radius option is invalid!');
end

% prepare labels
atom.name = obj.unit_cell.label(atom.idx);
% plot only the first word of every label
for ii = 1:nAtom
    labelTemp = strword(atom.name{ii},[1 2],true);
    label1 = labelTemp{1};
    %label2 = labelTemp{2};
    atom.text{ii}  = [label1 '(' num2str(atom.idx(ii)) ')_' num2str(ii)];
    %atom.ttip{ii}  = [label2 ' atom (' label1 ')' char(10) 'Unit cell:' char(10)];
end

% save atom coordinates into data
posDat = mat2cell([floor(pos);a2Idx],4,ones(1,size(pos,2)));


% shift positions
pos = bsxfun(@plus,pos,BV\param.shift);

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
    % number of atoms on the legend
    nLegend = sum(lSelect);
    % append color
    if strcmp(param.color,'auto')
        lDat.color = [lDat.color double(obj.unit_cell.color(:,lSelect))/255];
    else
        lDat.color = [lDat.color repmat(color/255,1,nLegend)];
    end
    % append type
    lDat.type = [lDat.type 3*ones(1,nLegend)];
    % append name
    lDat.name = [lDat.name repmat({'atom'},1,nLegend)];
    % append text
    lDat.text = [lDat.text obj.unit_cell.label(lSelect)];
    
    setappdata(hFigure,'legend',lDat);
    swplot.legend('on',hFigure);
end

% plot the atoms, text generated automatically
swplot.plot('type','ellipsoid','name','atom','position',pos,'R',radius,...
    'figure',hFigure,'color',color,'text','','legend',false,'label',lLabel,...
    'nmesh',param.nmesh,'tooltip',false,'data',posDat,'replace',param.replace,...
    'translate',param.translate,'zoom',param.zoom);

% save range
setappdata(hFigure,'range',struct('range',range,'unit',param.unit));

if param.tooltip
    swplot.tooltip('on',hFigure);
end

end