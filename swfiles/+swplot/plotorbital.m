function varargout = plotorbital(varargin)
% plots electron orbitals of Hydrogen
%
% SWPLOT.PLOTORBITAL('option1', value1, ...)
%
% hFigure = SWPLOT.PLOTORBITAL(...)
%
% The function plots selected electron orbitals of Hydrogen onto an swplot
% figure. The orbitals will be plotted on a selected central atom and a
% coordinate systems is defined by two positive ligand atoms (x-axis points
% towards the first ligand, y-axis points towards the second ligand).
% Orbitals will be drawn also on symmetry equivalent atoms by applying the
% the point group operations according to the crystal symmetry.
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
% mode      String, defines which orbital is plotted. Allowed strings and
%           the corresponding quantum numbers:
%
%               mode                n l m pm
%               's'         qNum = [1 0 0  0]
%               'p_x'       qNum = [2 1 1  1]
%               'p_y'       qNum = [2 1 0 -1]
%               'p_z'       qNum = [2 1 0  0]
%               'd_xy'      qNum = [3 2 2 -1]
%               'd_xz'      qNum = [3 2 1  1]
%               'd_yz'      qNum = [3 2 1 -1]
%               'd_z2'      qNum = [3 2 0  0]
%               'd_x2-y2'   qNum = [3 2 2  1]
%
%           where the quantum numbers are in a vector (n, l, m, {pm}),
%           where pm defines the optional linear combination of the +m and
%           -m orbitals (PSI is the wave function):
%               PSI = PSI(n,l,m) + pm*PSI(n,l,-m)
%           If pm is +1,-1 or m=0 the wave fuction is real, otherwise
%           complex.
% center    Center atom, defined by either a position vector (row vector
%           with coordinates in lattice units) or the label of the atom.
% ligand1   First positive ligand atom, defined by either a position vector
%           (row vector with coordinates in lattice units) or the label of
%           the atom.
% ligand2   Second positive ligand atom, defined by either a position
%           vector (row vector with coordinates in lattice units) or the
%           label of the atom.
% scale     Scaling factor for the size of the orbital relative to the
%           Hydrogen orbitals. Default value is 1.
% alpha     Transparency (alpha value) of the orbital, default value is 1.
% figure    Handle of the swplot figure. Default is the selected figure.
% color     Color of the orbital:
%               'auto'      All orbital gets the color of the ion.
%               'colorname' All orbital will have the same given color.
%               [R G B]     RGB code of the color that fix the color of all
%                           orbitals.
% color2    Second color for the orbital to differentiate between the lobes
%           of +/- PSI (wave function). Default is identical to the color
%           option.
% nPatch    Number of points in the surface mesh along the three
%           dimensions, default value is stored in
%           swpref.getpref('npatch').
% tooltip   If true, the tooltips will be shown when clicking on orbitals.
%           Default is true.
% shift     Column vector with 3 elements, all orbital positions will be
%           shifted by the given value. Default value is [0;0;0].
% replace   Replace previous orbital plot if true. Default is true.
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
% The name of the objects that are created called 'orbital'. To find the
% handles and the stored data on these objects, use e.g.
%
%   sObject = swplot.findobj(hFigure,'name','orbital')
%

% default values
%fontSize0 = swpref.getpref('fontsize',[]);
nMesh0    = swpref.getpref('nmesh',[]);
nPatch0   = swpref.getpref('npatch',[]);

inpForm.fname  = {'range' 'legend' 'label' 'scale' 'alpha' 'center'};
inpForm.defval = {[]      true     true    1       0.3     [0 0 0] };
inpForm.size   = {[-1 -2] [1 1]    [1 1]   [1 1]   [1 1]   [1 -8]  };
inpForm.soft   = {true    false    false   false   false   false   };

inpForm.fname  = [inpForm.fname  {'mode'  'color' 'nmesh' 'npatch' 'color2'}];
inpForm.defval = [inpForm.defval {'aniso' 'auto'  nMesh0  nPatch0  'auto'  }];
inpForm.size   = [inpForm.size   {[1 -4]  [1 -5]  [1 1]   [1 1]    [1 -7]  }];
inpForm.soft   = [inpForm.soft   {true    false   false   false    false   }];

inpForm.fname  = [inpForm.fname  {'figure' 'obj' 'unit' 'tooltip' 'copy'}];
inpForm.defval = [inpForm.defval {[]       []    'lu'   true      false }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1] [1 -6] [1 1]     [1 1] }];
inpForm.soft   = [inpForm.soft   {true     true  false  false     false }];

inpForm.fname  = [inpForm.fname  {'shift' 'replace' 'translate' 'zoom'}];
inpForm.defval = [inpForm.defval {[0;0;0] true      false       false }];
inpForm.size   = [inpForm.size   {[3 1]   [1 1]     [1 1]       [1 1] }];
inpForm.soft   = [inpForm.soft   {false   false     false       false }];

inpForm.fname  = [inpForm.fname  {'ligand1' 'ligand2'}];
inpForm.defval = [inpForm.defval {[0 0 0]   [0 0 0]  }];
inpForm.size   = [inpForm.size   {[1 -8]    [1 -9]   }];
inpForm.soft   = [inpForm.soft   {false     false    }];

param = sw_readparam(inpForm, varargin{:});

% find swplot figure
if isempty(param.figure)
    if isempty(param.obj)
        try
            hFigure  = swplot.activefigure;
        catch msg
            warning(msg.message)
            error('plotorbital:NoObj','The figure does not contain a SpinW object, use spinw.plot first!')
        end
    else
        hFigure  = swplot.activefigure('plot');
    end
else
    hFigure = param.figure;
end

% takes care of spinw object saved/loaded in/from figure
if isempty(param.obj)
    if isappdata(hFigure,'obj')
        obj = getappdata(hFigure,'obj');
    else
        error('plotorbital:NoObj','The figure does not contain a SpinW object, use spinw.plot first!')
    end
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
set(hFigure,'Name', 'SpinW: Electron orbitals');

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
    error('plotorbital:WrongInput','The given plotting range is invalid!');
end

switch param.unit
    case 'lu'
        rangelu = [floor(range(:,1)) ceil(range(:,2))];
    case 'xyz'
        % corners of the box
        cIdx    = cell2mat(arrayfun(@(M)bitget(M,1:3)+1,0:7,'UniformOutput',0)')';
        corners = BV\[range(1,cIdx(1,:));range(2,cIdx(2,:));range(3,cIdx(3,:))];
        rangelu = [min(corners,[],2) max(corners,[],2)];
        rangelu = [floor(rangelu(:,1)) ceil(rangelu(:,2))];
    otherwise
        error('plotorbital:WrongInput','The given unit string is invalid!');
end

% center


% atom data
mAtom  = obj.matom;

% generate bonds, but don't sort the bonds on DM interactions
[~, SI] = intmatrix(obj,'plotmode',true,'extend',false,'sortDM',false,'zeroC',false,'nExt',[1 1 1]);

% take care of return values
if nargout > 0
    varargout{1} = hFigure;
end

switch param.mode
    case 'aniso'
        if ~any(SI.aniso(:))
            % do nothing, just remove old plot
            % TODO
            return
        else
            mat = SI.aniso;
        end
    case 'g'
        mat = SI.g;
    otherwise
        error('plotorbital:WrongInput','The given mode string is invalid!');
end

% generate the  positions of the magnetic atoms
nMAtom = size(mAtom.r,2);

% generate positions from rangelu the inclusive range
pos = cell(1,3);
[pos{:}] = ndgrid(rangelu(1,1):rangelu(1,2),rangelu(2,1):rangelu(2,2),rangelu(3,1):rangelu(3,2));
pos = bsxfun(@plus,reshape(cat(4,pos{:}),[],3)',permute(mAtom.r,[1 3 2]));

% number of unit cells
nCell = size(pos,2);

% keep track of types of atoms
%aIdx = repmat(mAtom.idx,[nCell 1]);
% keep track of aniso matrix index
mIdx = repmat(1:nMAtom,[nCell 1]);
pos  = reshape(pos,3,[]);

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

pos  = pos(:,pIdx);
%aIdx = aIdx(pIdx);
mIdx = mIdx(pIdx);
%aIdx = aIdx(:)';
mIdx = mIdx(:)';
mat  = mat(:,:,mIdx);

% number of ellipse to plot
nEllipse = size(mat,3);

% color
if strcmp(param.color,'auto')
    color = double(obj.unit_cell.color(:,mAtom.idx(mIdx)));
else
    color = swplot.color(param.color);
end

if size(color,2) == 1
    color = repmat(color,1,numel(mIdx));
end

% save atom coordinates into data
%posDat = mat2cell(pos,3,ones(1,size(pos,2)));

% shift positions
pos = bsxfun(@plus,pos,BV\param.shift);

% length of shortest bond for scaling
if ~isempty(obj.coupling.atom1)
    bondlu = mAtom.r(:,obj.coupling.atom2(1))+double(obj.coupling.dl(:,1))-mAtom.r(:,obj.coupling.atom1(1));
    % minimum length of bond in Angstrom
    lxyz = min(sqrt(sum((BV*bondlu).^2,1)));
else
    lxyz = 3;
end


% plot ellipsoid
% remove zero ellipsoids
rmMat = permute(sumn(abs(mat),[1 2])==0,[1 3 2]);

mat   = mat(:,:,~rmMat);
pos   = pos(:,~rmMat);
color = color(:,~rmMat);
% calculating the main radiuses of the ellipsoid.
[V, Rell] = eigorth(mat,1e-5);
% creating positive definite matrix by adding constant to all
% eigenvalues.
maxR  = sqrt(max(sum(Rell.^2,1)));
switch param.mode
    case 'aniso'
        % large value --> short radius
        Rell = bsxfun(@minus,max(Rell,[],1),Rell);
    case 'g'
        
end
%Rell = ((Rell+param.radius1)/(maxR+param.radius1))*param.scale*min(lxyz);
Rell = (Rell/maxR)*param.scale*min(lxyz)+param.radius1;
% V*diag(R) vectorized
%V = bsxfun(@times,V,permute(Rell,[3 1 2]));
V = mmat(bsxfun(@times,V,permute(Rell,[3 1 2])),permute(V,[2 1 3]));

swplot.plot('type','ellipsoid','name','ion','position',pos,'text','',...
    'figure',hFigure,'legend','','color',color,'T',V,...
    'tooltip',false,'replace',param.replace,'nmesh',param.nmesh,...
    'data',{},'label',{},'translate',param.translate,...
    'zoom',param.zoom,'alpha',param.alpha);

if param.linewidth > 0
    % draw circles
    phi = linspace(0,2*pi,param.npatch+1);
    circle = [sin(phi);cos(phi);phi*0];
    % xy, xz and yz plane circles
    circle = [circle circle([1 3 2],:) circle([3 1 2],:)];
    
    %posc = permute(sum(bsxfun(@times,V,permute(circle,[3 1 4 2])),2),[1 3 4 2]);
    posc = permute(sum(bsxfun(@times,V,permute(circle,[3 1 4 2])),2),[1 3 4 2]);
    
    % convert to lu
    posc = reshape(BV\reshape(posc,3,[]),3,nEllipse,[],3);
    posc = reshape(permute(posc,[1 2 4 3]),3,3*nEllipse,[]);
    % shift to the atomic positions
    % color
    if strcmp(param.color2,'auto')
        color2 = color;
    else
        color2 = param.color2;
    end
    
    if size(color2,2)>1 && ~ischar(color2)
        color2 = repmat(color2,1,3);
    end
    
    posc = bsxfun(@plus,posc,repmat(pos,[1 3]));
    swplot.plot('type','line','name','ion_edge','position',posc,'text','',...
        'figure',hFigure,'legend','','color',color2,...
        'tooltip',false,'replace',param.replace,...
        'data',{},'label',{},'translate',param.translate,...
        'zoom',param.zoom,'alpha',param.alpha);
end

if param.tooltip
    swplot.tooltip('on',hFigure);
end

setappdata(hFigure,'range',struct('range',param.range,'unit',param.unit));

end