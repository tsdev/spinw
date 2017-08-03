function varargout = plotion(varargin)
% plots magnetic ion properties
%
% SWPLOT.PLOTION('option1', value1, ...)
%
% hFigure = SWPLOT.PLOTION(...)
%
% The function plots selected properties of magnetic ions stored in a SpinW
% object onto an swplot figure.
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
% mode      String, defines how the bond is plotted
%               'aniso'     Ellipsoid is plotted for single ion anisotropy.
%               'g'     	Ellipsoid is drawn for g-tensor.
% scale     Scaling factor for the size of the ellipsoid relative to the 
%           shortest bond length. Default value is 1/3.
% alpha     Transparency (alpha value) of the ellipsoid, default value is 
%           0.3.
% radius1   Minimum radius of the ellipsoid, default value is 0.08.
% lineWidth Line width in pt of the main circles surrounding the ellipsoid, 
%           if zero no circles are drawn. Default is 0.5.
% figure    Handle of the swplot figure. Default is the selected figure.
% legend    Whether to add the plot to the legend, default is true.
% color     Color of the ellipsoid:
%               'auto'      All ellipsoid gets the color of the ion.
%               'colorname' All ellipsoid will have the same given color.
%               [R G B]     RGB code of the color that fix the color of all
%                           ellipsoid.
% color2    Color of the main circles, default is 'auto' when the ellipses
%           will have the same color as the ellipsoids. Can be either a row
%           vector of RGB code or string of a color name.
% nMesh     Resolution of the ellipse surface mesh. Integer number that is
%           used to generate an icosahedron mesh with #mesh number of
%           additional triangulation, default value is stored in
%           swpref.getpref('nmesh')
% nPatch    Number of points on the curve for the arrows, default
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
% The name of the objects that are created called 'bond'. To find the
% handles and the stored data on these objects, use e.g.
%
%   sObject = swplot.findobj(hFigure,'name','bond')
%

% default values
%fontSize0 = swpref.getpref('fontsize',[]);
nMesh0    = swpref.getpref('nmesh',[]);
nPatch0   = swpref.getpref('npatch',[]);

inpForm.fname  = {'range' 'legend' 'label' 'scale' 'linewidth' 'alpha'};
inpForm.defval = {[]      true     true    1/3     0.5         0.3    };
inpForm.size   = {[-1 -2] [1 1]    [1 1]   [1 1]   [1 3]       [1 1]  };
inpForm.soft   = {true    false    false   false   false       false  };

inpForm.fname  = [inpForm.fname  {'mode'  'color' 'nmesh' 'npatch' 'color2'}];
inpForm.defval = [inpForm.defval {'aniso' 'auto'  nMesh0  nPatch0  'auto'  }];
inpForm.size   = [inpForm.size   {[1 -4]  [1 -5]  [1 1]   [1 1]    [1 -7]  }];
inpForm.soft   = [inpForm.soft   {true    false   false   false    false   }];

inpForm.fname  = [inpForm.fname  {'figure' 'obj' 'unit' 'tooltip' 'copy'}];
inpForm.defval = [inpForm.defval {[]       []    'lu'   true      false }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1] [1 -6] [1 1]     [1 1] }];
inpForm.soft   = [inpForm.soft   {true     true  false  false     false }];

inpForm.fname  = [inpForm.fname  {'shift' 'replace' 'radius1' 'translate' 'zoom'}];
inpForm.defval = [inpForm.defval {[0;0;0] true      0.08      false       false }];
inpForm.size   = [inpForm.size   {[3 1]   [1 1]     [1 1]     [1 1]       [1 1] }];
inpForm.soft   = [inpForm.soft   {false   false     false     false       false }];

param = sw_readparam(inpForm, varargin{:});

if isempty(param.figure)
    hFigure  = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

% take care of return values
if nargout > 0
    varargout{1} = hFigure;
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
set(hFigure,'Name', 'SpinW: Single ion properties');

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
    error('plotion:WrongInput','The given plotting range is invalid!');
end

switch param.unit
    case 'lu'
        rangelu = [floor(range(:,1)) ceil(range(:,2))];
    case 'xyz'
        % corners of the box
        cIdx    = cell2mat(arrayfun(@(M)bitget(M,1:3)+1,0:7,'UniformOutput',0)')';
        %corners = BV\[range(1,[1 2 2 1 1 2 2 1]);range(2,[1 1 2 2 1 1 2 2]);range(3,[1 1 1 1 2 2 2 2])];
        corners = BV\[range(1,cIdx(1,:));range(2,cIdx(2,:));range(3,cIdx(3,:))];
        rangelu = [min(corners,[],2) max(corners,[],2)];
        rangelu = [floor(rangelu(:,1)) ceil(rangelu(:,2))];
    otherwise
        error('plotion:WrongInput','The given unit string is invalid!');
end

% atom data
mAtom  = obj.matom;

% generate bonds, but don't sort the bonds on DM interactions
[~, SI] = intmatrix(obj,'plotmode',true,'extend',false,'sortDM',false,'zeroC',false,'nExt',[1 1 1]);

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
        error('plotion:WrongInput','The given mode string is invalid!');
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

% number of ellipse to plot
nEllipse = size(mat,3);

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