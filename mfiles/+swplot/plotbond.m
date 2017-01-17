function varargout = plotbond(varargin)
% plots magnetic bonds
%
% SWPLOT.PLOTBOND('option1', value1, ...)
%
% hFigure = SWPLOT.PLOTBOND(...)
%
% The function plots the magnetic bonds stored in a SpinW object onto an
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
% mode      String, defines how the bond is plotted
%               'cylinder'  Bonds are plotted as cylinders (default).
%               'arrow'     Bonds are plotted as arrows (default if DM
%                           interactions are non-zero).
%               'line'      Bonds are plotted as continuous/dashed lines
%                           depending on the label of the corresponding
%                           matrix (dashed line is used if the matrix
%                           label ends with '-', otherwise continuous).
%               '--'        Bonds are plotted as dashed lines.
%               '-'         Bonds are plotted as lines.
% mode2     String, defines what is plotted on the bond:
%               'none'      Don't plot anything on the bond (default).
%               'antisym'   Plot the antisymmetric part (DM vector) of the 
%                           exchange at the middle point of the bond
%                           (default if DM vectors are non-zero).
%               'sym'       Plot the symmetric exchange at the middle
%                           of the bond as an ellipsoid.
% linewidth Defines the bond radius if it is drawn by a line:
%               'fix'       All line will have a width given by linewidth0.
%                           Default value.
%               'auto'      Lines will have a width that is depending on
%                           the exchange matrix on the bond: W~sum(abs(J)),
%                           where the largest line width on the strongest
%                           bond is given by linewidth0.
% linewidth0 Line width in pt used to draw the bond if 'mode' is 'line',
%           '--' or '-'. Default value is 0.5.
% zero      If true, bonds with zero exchange matrix will be plotted as
%           well. Default is true.
% radius0   Radius of the cylinder, default value is 0.05.
% radius2   Constant atom radius, default value is 0.3 Angstrom.
% radius    Defines the atom radius (important for arrow bonds, to avoid
%           overlap with the spheres of the atoms):
%               'fix'       Sets the radius of all atoms to the value
%                           given by radius2.
%               'auto'      use radius data from database based on the atom
%                           label multiplied by radius2 value.
% ang       Angle of the arrow head in degree units, default is 30 degree.
% lHead     Length of the arrow head, default value is 0.5.
% scale     Scaling factor for the lenght of the DM vector or the size of
%           the ellipsoid. Default value is 1.
% figure    Handle of the swplot figure. Default is the selected figure.
% legend    Whether to add the plot to the legend, default is true.
% color     Color of the bonds:
%               'auto'      All bonds get the stored color.
%               'colorname' All bonds will have the same given color.
%               [R G B]     RGB code of the color that fix the color of all
%                           bonds.
% color2    Color of the ellipse or DM vector on the bond:
%               'auto'      All object get the color of the bond.
%               'colorname' All object will have the same given color.
%               [R G B]     RGB code of the color that fix the color of all
%                           object.
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
range0    = [0 1;0 1;0 1];

inpForm.fname  = {'range' 'legend' 'label' 'zero' 'scale' 'radius0' 'mode2' 'linewidth' };
inpForm.defval = {range0  true     true    true   1       0.05      []      'fix'       };
inpForm.size   = {[-1 -2] [1 1]    [1 1]   [1 1]  [1 1]   [1 1]     [1 -7]  [1 -9]      };
inpForm.soft   = {false   false    false   false  false   false     true    false       };

inpForm.fname  = [inpForm.fname  {'radius' 'mode' 'color' 'nmesh' 'npatch' 'linewidth0' }];
inpForm.defval = [inpForm.defval {'auto'   []     'auto'  nMesh0  nPatch0  0.5          }];
inpForm.size   = [inpForm.size   {[1 -3]   [1 -4] [1 -5]  [1 1]   [1 1]    [1 1]        }];
inpForm.soft   = [inpForm.soft   {false    true  false   false   false    false         }];

inpForm.fname  = [inpForm.fname  {'figure' 'obj' 'rangeunit' 'tooltip' 'radius'}];
inpForm.defval = [inpForm.defval {[]       []    'lu'        true      'auto'  }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1] [1 -6]      [1 1]     [1 -8]  }];
inpForm.soft   = [inpForm.soft   {true     true  false       false     false   }];

inpForm.fname  = [inpForm.fname  {'shift' 'replace' 'arrow' 'ang' 'lHead' 'radius2'}];
inpForm.defval = [inpForm.defval {[0;0;0] true      []      30    0.5     0.3      }];
inpForm.size   = [inpForm.size   {[3 1]   [1 1]     [1 1]   [1 1] [1 1]   [1 1]    }];
inpForm.soft   = [inpForm.soft   {false   false     true    false false   false    }];

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
set(hFigure,'Name', 'SpinW: Magnetic bonds');

% change range, if the number of unit cells are given
if numel(param.range) == 3
    param.range = [ zeros(3,1) param.range(:)];
elseif numel(param.range) ~=6
    error('plotbond:WrongInput','The given plotting range is invalid!');
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
        error('plotbond:WrongInput','The given rangeunit string is invalid!');
end

% atom data
mAtom  = obj.matom;

% generate bonds, but don't sort the bonds on DM interactions
[SS, ~] = intmatrix(obj,'plotmode',true,'extend',false,'sortDM',false,'zeroC',param.zero,'nExt',[1 1 1]);

if isempty(SS.all)
    warning('plotbond:EmptyPlot','No bonds to plot!')
end

coupling        = struct;
coupling.dl     = double(SS.all(1:3,:));
coupling.atom1  = SS.all(4,:);
coupling.atom2  = SS.all(5,:);
coupling.matidx = SS.all(15,:);
coupling.idx    = SS.all(16,:);

% matrix values
mat   = reshape(SS.all(6:14,:),3,3,[]);
% DM interaction
matDM = (mat-permute(mat,[2 1 3]))/2;
% keep symmetric part of the interactions
matSym = (mat+permute(mat,[2 1 3]))/2;
matDM  = permute(cat(2,matDM(2,3,:),matDM(3,1,:),matDM(1,2,:)),[2 3 1]);

% are there non-zero DM vectors
isDM = any(matDM(:));

% select bond type to plot
if isempty(param.mode)
    if isDM
        param.mode = 'arrow';
    else
        param.mode = 'cylinder';
    end
end

if isempty(param.mode2)
    if isDM
        param.mode2 = 'antisym';
    else
        param.mode2 = 'none';
    end
end
  

% generate the  positions of the bonds
% generate positions from rangelu the inclusive range
pos = cell(1,3);
[pos{:}] = ndgrid(rangelu(1,1):rangelu(1,2),rangelu(2,1):rangelu(2,2),rangelu(3,1):rangelu(3,2));
pos  = reshape(cat(4,pos{:}),[],3)';
pos1 = bsxfun(@plus,pos,permute(mAtom.r(:,coupling.atom1),[1 3 2]));
pos2 = bsxfun(@plus,pos,permute(mAtom.r(:,coupling.atom2)+coupling.dl,[1 3 2]));

% number of unit cells
nCell = size(pos,2);

% generate matrix indices
matidx = repmat(coupling.matidx,[nCell 1]);

% generate matrix values
matSym = reshape(repmat(permute(matSym,[1 2 4 3]),[1 1 nCell 1]),3,3,[]);
matDM  = reshape(repmat(permute(matDM,[1 3 2]),[1 nCell 1]),3,[]);

pos1  = reshape(pos1,3,[]);
pos2  = reshape(pos2,3,[]);

% cut out the bonds that are out of range
switch param.rangeunit
    case 'lu'
        % lower range<=L<= upper range
        pIdx1 = all(bsxfun(@ge,pos1,range(:,1)) & bsxfun(@le,pos1,range(:,2)),1);
        pIdx2 = all(bsxfun(@ge,pos2,range(:,1)) & bsxfun(@le,pos2,range(:,2)),1);
        pIdx  = all([pIdx1;pIdx2],1);
    case 'xyz'
        % convert to xyz
        posxyz1 = BV*pos1;
        posxyz2 = BV*pos2;
        pIdx1   = all(bsxfun(@ge,posxyz1,range(:,1)) & bsxfun(@le,posxyz1,range(:,2)),1);
        pIdx2   = all(bsxfun(@ge,posxyz2,range(:,1)) & bsxfun(@le,posxyz2,range(:,2)),1);
        pIdx    = all([pIdx1;pIdx2],1);
end

if ~any(pIdx)
    warning('plotbond:EmptyPlot','There are no bonds in the plotting range!')
    return
end

pos1   = pos1(:,pIdx);
pos2   = pos2(:,pIdx);
matidx = matidx(pIdx);
matSym = matSym(:,:,pIdx);
matDM  = matDM(:,pIdx);

% shift positions
pos1 = bsxfun(@plus,pos1,param.shift);
pos2 = bsxfun(@plus,pos2,param.shift);

% color
if strcmp(param.color,'auto')
    color = double(obj.matrix.color(:,matidx));
else
    color = swplot.color(param.color);
end

% save original matrix values into data
mat0   = obj.matrix.mat(:,:,matidx);
matDat = mat2cell(mat0,3,3,ones(1,numel(matidx)));

% legend label
lLabel = obj.matrix.label(matidx);

lineStyle0 = '-';
switch param.mode
    case 'cylinder'
        type0 = 'cylinder';
    case 'arrow'
        type0 = 'arrow';
        % shorten the bonds
        % radius
        switch param.radius
            case 'auto'
                radius = sw_atomdata(obj.unit_cell.label,'radius');
                aidx1 = repmat(mAtom.idx(1,coupling.atom1),[nCell 1]);
                aidx1 = aidx1(pIdx);
                aidx2 = repmat(mAtom.idx(1,coupling.atom2),[nCell 1]);
                aidx2 = aidx2(pIdx);
                rad1 = radius(aidx1)*param.radius2;
                rad2 = radius(aidx2)*param.radius2;
            case 'fix'
                rad1 = param.radius2;
                rad2 = param.radius2;
            otherwise
                error('plotmag:WrongInput','The given radius string is invalid!');
        end
        % shorten
        dpos = pos2-pos1;
        % length of bond
        dposxyz = BV*dpos;
        lxyz = sqrt(sum(dposxyz.^2,1));
        
        pos1s = pos1 + bsxfun(@times,rad1./lxyz,dpos);
        pos2s = pos2 - bsxfun(@times,rad2./lxyz,dpos);
    case 'line'
        type0 = 'line';
        % change linestyle based on the matrix label
        % 1 for continuous line
        % 2 for dashed line
        lineStyle0 = (cellfun(@(C)C(end),lLabel)=='-')+1;
    case '--'
        type0 = 'line';
        lineStyle0 = '--';
    case '-'
        type0 = 'line';
    otherwise
        error('plotbond:WrongInput','The given mode string is illegal!')
end
    
switch param.linewidth
    case 'fix'
        lineWidth = param.linewidth0;
    case 'auto'
        absmat = permute(sumn(abs(mat),[1 2]),[1 3 2]);
        lineWidth = absmat/max(absmat)*param.linewidth0;
    otherwise
        error('plotbond:WrongInput','The given lineStyle string is illegal!')
end

% plot bond vectors
swplot.plot('type',type0,'name','bond','position',cat(3,pos1s,pos2s),'text','',...
    'figure',hFigure,'legend',lLabel,'color',color,'R',param.radius0,...
    'tooltip',false,'replace',param.replace,'lineStyle',lineStyle0,...
    'data',matDat,'label',{},'nmesh',param.nmesh,'ang',param.ang,...
    'lineWidth',lineWidth,'lHead',param.lHead);

% plot ellipse/arrow on top of bond
switch param.mode2
    case 'none'
        % do nothing
    case 'antisym'
        % plot arrows
        % TODO
        swplot.plot('type','arrow','name','bond_DM','position',cat(3,(pos1+pos2)/2,(pos1+pos2)/2+1),'text','',...
            'figure',hFigure,'legend',lLabel,'color',color,'R',param.radius0,...
            'tooltip',false,'replace',param.replace,'lineStyle',lineStyle0,...
            'data',matDat,'label',{},'nmesh',param.nmesh,'ang',param.ang,...
            'lineWidth',lineWidth,'lHead',param.lHead);
        
    case 'sym'
    otherwise
        error('plotbond:WrongInput','The given mode2 string is invalid!');
end

% legend data
lDat = getappdata(hFigure,'legend');

if param.replace
    % remove old legend entries
    lIdx = ~ismember(lDat.name,'bond');
    lDat.color = lDat.color(:,lIdx);
    lDat.type  = lDat.type(:,lIdx);
    lDat.name  = lDat.name(:,lIdx);
    lDat.text  = lDat.text(:,lIdx);
    % redraw legend
    swplot.legend('refresh',hFigure);
end

if param.legend
    % append color
    lDat.color = [lDat.color double(obj.matrix.color)/255];
    % append type
    lDat.type = [lDat.type (cellfun(@(C)C(end),obj.matrix.label) == '-')+1];
    % append name
    lDat.name = [lDat.name repmat({'bond'},1,numel(obj.matrix.label))];
    % append text
    lDat.text = [lDat.text obj.matrix.label];
    
    setappdata(hFigure,'legend',lDat);
    swplot.legend('on',hFigure);
end

if nargout > 0
    varargout{1} = hFigure;
end

if param.tooltip
    swplot.tooltip('on',hFigure);
end

end