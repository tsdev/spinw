function varargout = plotorbital(varargin)
% plots electron orbitals of hydrogen
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
% mode      String, determines the kind of orbital to plot. The only
%           allowed value is 'hydrogen' for plotting hydrogen atomic
%           orbitals.
% type      String, defines which orbital is plotted. Allowed strings and
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
%           hydrogen orbitals. Default value is 0.15.
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
inpForm.defval = {[]      true     true    0.15    1       [0 0 0] };
inpForm.size   = {[-1 -2] [1 1]    [1 1]   [1 1]   [1 1]   [1 -8]  };
inpForm.soft   = {true    false    false   false   false   false   };

inpForm.fname  = [inpForm.fname  {'mode'     'color' 'nmesh' 'npatch' 'color2'}];
inpForm.defval = [inpForm.defval {'hydrogen' 'auto'  nMesh0  nPatch0  'auto'  }];
inpForm.size   = [inpForm.size   {[1 -4]     [1 -5]  [1 1]   [1 1]    [1 -7]  }];
inpForm.soft   = [inpForm.soft   {true       false   false   false    false   }];

inpForm.fname  = [inpForm.fname  {'figure' 'obj' 'unit' 'tooltip' 'copy'}];
inpForm.defval = [inpForm.defval {[]       []    'lu'   true      false }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1] [1 -6] [1 1]     [1 1] }];
inpForm.soft   = [inpForm.soft   {true     true  false  false     false }];

inpForm.fname  = [inpForm.fname  {'shift' 'replace' 'translate' 'zoom'}];
inpForm.defval = [inpForm.defval {[0;0;0] true      false       false }];
inpForm.size   = [inpForm.size   {[3 1]   [1 1]     [1 1]       [1 1] }];
inpForm.soft   = [inpForm.soft   {false   false     false       false }];

inpForm.fname  = [inpForm.fname  {'ligand1' 'ligand2' 'type' }];
inpForm.defval = [inpForm.defval {[1 0 0]   [0 1 0]   'd_yz' }];
inpForm.size   = [inpForm.size   {[1 -9]    [1 -10]   [1 -11]}];
inpForm.soft   = [inpForm.soft   {false     false     false  }];

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

% is the position vector for the xy-axis are given directly
posvec = false;

% positions from string
pPos0 = {param.center param.ligand1 param.ligand2};
pStr0 = {'center' 'ligand1' 'ligand2'};
dat   = struct('pos',zeros(3),'idx',0);
for ii = 1:3
    pPos = pPos0{ii};
    pStr = pStr0{ii};
    if ischar(pPos)
        % a symmetry position is defined with an atom label string
        atomIdx = cellfun(@(C)~isempty(C),strfind(obj.unit_cell.label,pPos));
        if sum(atomIdx)>1
            error('plotorbital:WrongLabel','The label of the %s atom is not unique!',pStr)
        elseif sum(atomIdx)==0
            error('plotorbital:WrongLabel','The given label for the %s atom does not exists!',pStr)
        end
        
        if ii==1
            % just save the first symmetry positon for center atom
            dat.idx = find(atomIdx);
            dat.pos(:,ii) = obj.unit_cell.r(:,atomIdx);
        else
            % find the closest ligand position from symmetry equivalent
            % positions
            rLigand = obj.atom.r(:,obj.atom.idx==find(atomIdx));
            if ii == 3
                % remove the 1st ligand position from the list
                rLigand(:,all(bsxfun(@eq,dat.pos(:,2),rLigand),1)) = [];
            end
            % find the first shortest distance
            dat.pos(:,ii) = rLigand(:,findmin(sum((BV*bsxfun(@minus,dat.pos(:,1),rLigand)).^2,1),1));
        end
        
    elseif numel(pPos) == 1
        % atom position for index in spinw.atom list
        dat.pos(:,ii) = obj.atom.r(:,pPos);
        if ii == 1
            dat.idx       = obj.atom.idx(pPos);
        end
    elseif numel(pPos) == 3
        if ii == 1
            % shift center position into the first cell
            dat.pos(:,ii) = mod(pPos',1);
        else
            % directly give the x-axis for the ligands
            dat.pos(:,ii) = dat.pos(:,1) + pPos';
            posvec = true;
        end
    else
        error('plotorbital:WrongInput','The given input for %s atom is wrong!',pStr)
    end
    
end

% find shortest distance by shifting the unit cells
% generate -1,0,+1 3D grid stored in R
% only when the position vectors are not given directly
if ~posvec
    R = cell(1,3);
    [R{:}] = ndgrid(-1:1,-1:1,-1:1);
    R = reshape(cat(4,R{:}),[],3)';
    % check: swplot.plot('type','ellipsoid','R',0.1,'position',permute(R,[2 1 3]))
    
    for ii = 2:3
        % only shift position if it is significantly shorter than the original
        % the given limit is 0.1 Angstrom
        R0 = sqrt(sum((BV*(dat.pos(:,1)-dat.pos(:,ii))).^2));
        posNew = dat.pos(:,ii)+R(:,findmin(sum((BV*bsxfun(@minus,dat.pos(:,1)-dat.pos(:,ii),R)).^2,1),1));
        R1 = sqrt(sum((BV*(dat.pos(:,1)-posNew)).^2));
        if R1 < R0 - 0.1
            dat.pos(:,ii) = posNew;
        end
    end
end

% determine the color of the orbitals
% first color
if strcmp(param.color,'auto')
    if dat.idx > 0
        % color the orbital according to the center atom, use single color
        color = double(obj.unit_cell.color(:,dat.idx));
    else
        % use default color: red
        color = [255;0;0];
    end
elseif ischar(param.color)
    color = swplot.color(param.color);
elseif numel(param.color)==3
    color = param.color(:);
else
    error('plotorbital:WrongInput','The format of the given color option is invalid!')
end

% second color
if strcmp(param.color2,'auto')
    % find complementer color
    color(:,2) = complementer(color(:,1));
elseif ischar(param.color2)
    color(:,2) = swplot.color(param.color2);
elseif numel(param.color2)==3
    color = param.color2(:);
else
    error('plotorbital:WrongInput','The format of the given color2 option is invalid!')
end

% determine the T0 Cartesian coordinate system for the first atom
posXYZ = BV*dat.pos;
dat.T0 = sw_cartesian([diff(posXYZ(:,1:2),1,2),diff(posXYZ(:,[1 3]),1,2)]);

% find the symmetry equivalent positions
[dat.allpos,~,symInfo] = swsym.position(obj.lattice.sym,dat.pos(:,1));
dat.R = mmat(symInfo.opmove,inv(symInfo.opmove(:,:,1)));
dat.R = mmat(BV,mmat(dat.R,inv(BV)));

% rotate the basis vectors on each symmetry equivalent position
dat.T = mmat(dat.R,dat.T0);

% test: draw vectors towards the ligands
% r1 = dat.pos(:,[1 1]);
% r2 = dat.pos(:,[2 3]);
% swplot.plot('type','arrow','position',cat(3,r1,r2),'color',{'red','blue'})


% take care of return values
if nargout > 0
    varargout{1} = hFigure;
end

switch param.mode
    case 'hydrogen'
        % do nothing yet
    otherwise
        error('plotorbital:WrongInput','The given mode string is invalid!');
end

% generate positions from rangelu the inclusive range
pos = cell(1,3);
[pos{:}] = ndgrid(rangelu(1,1):rangelu(1,2),rangelu(2,1):rangelu(2,2),rangelu(3,1):rangelu(3,2));
pos = bsxfun(@plus,reshape(cat(4,pos{:}),[],3)',permute(dat.allpos,[1 3 2]));

% number of unit cells
nCell = size(pos,2);

% keep track of index of orbital in the original cell
oIdx = repmat(1:size(dat.allpos,2),[nCell 1]);
pos  = reshape(pos,3,[]);
oIdx = oIdx(:);

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
    warning('plotorbital:EmptyPlot','There are no orbitals in the plotting range!')
    return
end

pos  = pos(:,pIdx);
oIdx = oIdx(pIdx);

% generate all transformation matrix
T = dat.T(:,:,oIdx);

% shift positions
pos = bsxfun(@plus,pos,BV\param.shift);

% draw the orbitals
swplot.plot('type','orbital','qLabel',param.type,'position',pos,'T',T,...
    'scale',param.scale,'nPatch',param.npatch,'color',permute(color,[1 3 2]),...
    'alpha',param.alpha);

if param.tooltip
    swplot.tooltip('on',hFigure);
end

setappdata(hFigure,'range',struct('range',param.range,'unit',param.unit));

end

function RGB = complementer(RGB)
% find complementer RGB color
%
% RGB = swplot.complementer(RGB)
%

HSV    = rgb2hsv(RGB(:)'/255);
HSV(1) = mod(HSV(1)+1/2,1);
RGB    = hsv2rgb(HSV)*255;

end