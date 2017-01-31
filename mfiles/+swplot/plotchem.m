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
%           0.3 for polyhedron.
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
%           center. Default is false.
% zoom      If true, figure will be automatically zoomed to the ideal size.
%           Default is false.
% copy      If true, a hardcopy of the spinw object will be sved in the
%           figure data, otherwise just the handle of the spinw object, 
%           thus the figure can be updated when the spin object changed.
%           Default value is false. 
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

inpForm.fname  = {'range' 'legend' 'label' 'unit' 'mode' 'atom1' 'atom2' 'copy'};
inpForm.defval = {range0  true     true    'lu'   'poly' 1       2       false };
inpForm.size   = {[-1 -2] [1 1]    [1 1]   [1 -3] [1 -4] [1 -5]  [1 -6]  [1 1] };
inpForm.soft   = {false   false    false   false  false  false   false   false };

inpForm.fname  = [inpForm.fname  {'limit' 'alpha' 'color' 'nmesh' 'npatch'}];
inpForm.defval = [inpForm.defval {6       []      'auto'  nMesh0  nPatch0 }];
inpForm.size   = [inpForm.size   {[1 -7]  [1 1]   [1 -8]  [1 1]   [1 1]   }];
inpForm.soft   = [inpForm.soft   {false   true    false   false   false   }];

inpForm.fname  = [inpForm.fname  {'figure' 'obj' 'color2' 'tooltip' 'radius0'}];
inpForm.defval = [inpForm.defval {[]       []    'auto'   true      0.03     }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1] [1 -9]  [1 1]      [1 1]    }];
inpForm.soft   = [inpForm.soft   {true     true  false   false      false    }];

inpForm.fname  = [inpForm.fname  {'shift' 'replace' 'translate' 'zoom' 'extend' }];
inpForm.defval = [inpForm.defval {[0;0;0] true      false        false  2       }];
inpForm.size   = [inpForm.size   {[3 1]   [1 1]     [1 1]        [1 1]  [1 1]   }];
inpForm.soft   = [inpForm.soft   {false   false     false        false  false   }];

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
    warning('plotchem:WrongInput','No SpinW object to plot!');
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

% convert atom label to index
if ischar(param.atom1)
    atom1 = find(cellfun(@isempty,strfind(obj.unit_cell.label,param.atom1))==false);
else
    atom1 = param.atom1;
end

if ischar(param.atom2)
    atom2 = find(cellfun(@isempty,strfind(obj.unit_cell.label,param.atom2))==false);
else
    atom2 = param.atom2;
end

if isempty(atom1) || isempty(atom2)
    warning('plotchem:EmptyPlot','There are no atoms in the plotting range!')
    return
end

% find the indices of atom1 and atom2 in the atom list
atom1idx = find(ismember(atom.idx,atom1));
atom2idx = find(ismember(atom.idx,atom2));


% select atom type to plot
switch param.mode
    case 'poly'
    case 'bond'
    otherwise
        error('plotchem:WrongInput','The given mode string is invalid!');
end

%nAtom = size(atom.r,2);

% generate atom1 positions from rangelu the inclusive range
pos1 = cell(1,3);
[pos1{:}] = ndgrid(rangelu(1,1):rangelu(1,2),rangelu(2,1):rangelu(2,2),rangelu(3,1):rangelu(3,2));
pos1      = bsxfun(@plus,reshape(cat(4,pos1{:}),[],3)',permute(atom.r(:,atom1idx),[1 3 2]));
nCell1    = size(pos1,2);
% keep track of types in atom1 list
atom1idx = repmat(atom1idx,[nCell1 1]);
pos1     = reshape(pos1,3,[]);

% generate atom1 positions from rangelu the inclusive range
pos2 = cell(1,3);
% extend range for second atom by 2 unit cell in all directions
rangeluext = [rangelu(:,1)-param.extend rangelu(:,2)+param.extend];
[pos2{:}]  = ndgrid(rangeluext(1,1):rangeluext(1,2),rangeluext(2,1):rangeluext(2,2),rangeluext(3,1):rangeluext(3,2));
pos2       = bsxfun(@plus,reshape(cat(4,pos2{:}),[],3)',permute(atom.r(:,atom2idx),[1 3 2]));
nCell2     = size(pos2,2);
% keep track of types in atom1 list
atom2idx  = repmat(atom2idx,[nCell2 1]);
pos2      = reshape(pos2,3,[]);

% cut out the atom1 positions that are out of range
switch param.unit
    case 'lu'
        % L>= lower range, L<= upper range
        pIdx = all(bsxfun(@ge,pos1,range(:,1)) & bsxfun(@le,pos1,range(:,2)),1);
    case 'xyz'
        % convert to xyz
        posxyz = BV*pos1;
        pIdx = all(bsxfun(@ge,posxyz,range(:,1)) & bsxfun(@le,posxyz,range(:,2)),1);
end

if ~any(pIdx)
    warning('plotchem:EmptyPlot','There are no atoms in the plotting range!')
    return
end

pos1     = pos1(:,pIdx);
atom1idx = atom1idx(pIdx);
atom1idx = atom1idx(:)';
atom2idx = atom2idx(:)';

% atom position matrices
pos1 = repmat(permute(pos1,[2 3 1]),[1 numel(atom2idx)]);
pos2 = repmat(permute(pos2,[3 2 1]),[numel(atom1idx) 1]);
% atom index matrices
atom1idx = repmat(atom1idx',[1 numel(atom2idx)]);
%atom2idx = repmat(atom2idx ,[size(atom1idx,1) 1]);
% calculate bond vectors in lu
dist = bsxfun(@minus,pos1,pos2);
% convert to xyz and calculate length
dist = sqrt(sum(sum(bsxfun(@times,permute(BV,[3 4 1 2]),dist),4).^2,3));
% sort bond length
[dist,sIdx] = sort(dist,2);
% identify unique atom1
sortidx1 = repmat((1:size(atom1idx,1))',[1 size(atom1idx,2)]);

% add extra line for fake atoms


% find bonds that fulfill the given limits
limit = param.limit;
if numel(limit) == 1
    % fixed number of bonds
    if (limit-round(limit))~=0
        error('plotchem:WrongInput','If limit option is a single number, it has to be integer!')
    end
    if limit>size(sIdx,2)
        error('plotchem:WrongInput','Increase the extend option (2) to find all furher neighbors!')
    end
    % convert into linear index from sort
    sortidx1 = sortidx1(:,1:limit);
    sIdx     = sIdx(:,1:limit);
    
    linIdx   = sub2ind(size(dist),sortidx1,sIdx);
    % select the atoms that fulfill the limits
    atom1idx = atom1idx(linIdx);
    %atom2idx = atom2idx(linIdx);
    % also select the position vectors
    pos1 = cat(3,pos1(linIdx),pos1(linIdx+numel(dist)),pos1(linIdx+2*numel(dist)));
    pos2 = cat(3,pos2(linIdx),pos2(linIdx+numel(dist)),pos2(linIdx+2*numel(dist)));

elseif numel(limit)==2
    nBond1 = sum(dist>=limit(1) & dist<=limit(2),2);
    
    if any(diff(nBond1))
        % polyhedra with different number of vertices
        % loop over all atoms
    else
        % all polyhedra have the same number of vertices
        sIdx = sIdx';
        sortidx1 = sortidx1';
        idx0 = (dist')>=limit(1) & (dist')<=limit(2);
        sIdx = reshape(sIdx(idx0),nBond1(1),[])';
        sortidx1 = reshape(sortidx1(idx0),nBond1(1),[])';
    end
    
    linIdx   = sub2ind(size(dist),sortidx1,sIdx);
    % select the atoms that fulfill the limits
    atom1idx = atom1idx(linIdx);
    %atom2idx = atom2idx(linIdx);
    % also select the position vectors
    pos1 = cat(3,pos1(linIdx),pos1(linIdx+numel(dist)),pos1(linIdx+2*numel(dist)));
    pos2 = cat(3,pos2(linIdx),pos2(linIdx+numel(dist)),pos2(linIdx+2*numel(dist)));
end

% shift positions
pos1 = bsxfun(@plus,pos1,permute(BV\param.shift,[2 3 1]));
pos2 = bsxfun(@plus,pos2,permute(BV\param.shift,[2 3 1]));

switch param.mode
    case {'poly' 'polyhedron'}
        if isempty(param.alpha)
            param.alpha = 0.3;
        end
        
        % reshape pos to [3 nObject nVertex]
        pos = permute(pos2,[3 1 2]);

        % color of the polyhedron faces
        if strcmp(param.color,'auto')
            color = double(obj.unit_cell.color(:,atom.idx(atom1idx(:,1))));
        else
            color = swplot.color(param.color);
        end
        
        % plot the atoms, text generated automatically
        swplot.plot('type','polyhedron','name','chem','position',pos,...
            'figure',hFigure,'color',color,'text','','legend',false,'label','',...
            'nmesh',param.nmesh,'tooltip',false,'data',{},'replace',param.replace,...
            'translate',param.translate,'zoom',param.zoom,'alpha',param.alpha);
        
        if strcmp(param.color2,'auto')
            % use the first color
            color2 = color(:,1);
            % change edge color of the last object (the poly we added)
            sObject = getappdata(hFigure,'objects');
            set(sObject(end).handle,'EdgeColor',color2/255);
        elseif strcmp(param.color2,'none')
            % do nothing
        else
            color2 = swplot.color(param.color2);
            % change edge color of the last object (the poly we added)
            sObject = getappdata(hFigure,'objects');
            set(sObject(end).handle,'EdgeColor',color2);
        end
        
    case 'bond'
        if isempty(param.alpha)
            param.alpha = 1;
        end

        % reshape pos to [3 nObject nVertex]
        pos = cat(3,reshape(permute(pos1,[3 1 2]),3,[]),reshape(permute(pos2,[3 1 2]),3,[]));

        if strcmp(param.color,'auto')
            color = double(obj.unit_cell.color(:,atom.idx(atom1idx(:))));
        else
            color = swplot.color(param.color);
        end
        
        % plot the atoms, text generated automatically
        swplot.plot('type','cylinder','name','chem','position',pos,...
            'figure',hFigure,'color',color,'text','','legend',false,'label','',...
            'nmesh',param.nmesh,'tooltip',false,'data',{},'replace',param.replace,...
            'translate',param.translate,'zoom',param.zoom);
    otherwise
        error('plotchem:WrongInput','The mode string is invalid!');
end
% 
% % legend label
% lLabel = repmat(atom.name,[nCell 1]);
% lLabel = lLabel(pIdx);
% 
% % legend data
% lDat = getappdata(hFigure,'legend');
% 
% if param.replace
%     % remove old legend entries
%     lIdx = ~ismember(lDat.name,'atom');
%     lDat.color = lDat.color(:,lIdx);
%     lDat.type  = lDat.type(:,lIdx);
%     lDat.name  = lDat.name(:,lIdx);
%     lDat.text  = lDat.text(:,lIdx);
%     setappdata(hFigure,'legend',lDat);
%     % redraw legend
%     swplot.legend('refresh',hFigure);
% end
% 
% if param.legend
%     % append color
%     lDat.color = [lDat.color double(obj.unit_cell.color)/255];
%     % append type
%     lDat.type = [lDat.type 3*ones(1,obj.natom)];
%     % append name
%     lDat.name = [lDat.name repmat({'atom'},1,obj.natom)];
%     % append text
%     lDat.text = [lDat.text obj.unit_cell.label];
%     
%     setappdata(hFigure,'legend',lDat);
%     swplot.legend('on',hFigure);
% end

% save range
setappdata(hFigure,'range',struct('range',param.range,'unit',param.unit));

if param.tooltip
    swplot.tooltip('on',hFigure);
end

end