function varargout = plotbond(varargin)
% plots bonds
% 
% ### Syntax
% 
% `swplot.plotbond(Name,Value)`
% 
% `hFigure = swplot.plotbond(Name,Value)`
%
% ### Description
% 
% `swplot.plotbond(Name,Value)` plots the magnetic bonds stored in
% [spinw.coupling]. It can plot bonds using different styles, such as
% arrows, lines or cylinders and allows controlling the the line style,
% cylinder radius, color, etc. The command can also plot cylinders with
% radii that are dependent on the strength of the exchange interaction
% assigned to the particular bond.
% 
% ### Name-Value Pair Arguments
% 
% `'mode'`
% : String that defines the style of the bond:
%   * `'cylinder'`   bonds are plotted as cylinders connecting the two atoms
%                   (default),
%   * `'arrow'`      bonds are plotted as arrows (default if DM
%                   interactions are present),
%   * `'line'`       bonds are plotted as lines,
%   * `'empty'`      no bonds are plotted.
% 
% `'mode2'`
% : String that defines what object is plotted on the middle of the bond:
%   * `'none'`      don't plot anything on the bond (default),
%   * `'antisym'`   plot the antisymmetric part (DM vector) of the 
%                   exchange on the middle point of the bond
%                   (default if DM vectors are non-zero),
%   * `'sym'`       plot an ellipsoid (corresponding to the symmetric
%                   exchange) on the middle of the bond.
% 
% `'sign'`
% : String that defines how the ellipsoids are generated for exchange
%   matrices that contain both negative and positive eigenvalues.
%   Possible values are:
%   * `'abs'`       The absolute value of the eigenvalues is used,
%                   this works nicely except that AFM and FM values
%                   will have the same radius (default).
%   * `'min'`       If there is a negative eigenvalue, it is
%                   shifted to zero with all other egeinvalues
%                   equally. This works nicely to emphasize AFM
%                   type values in the exchange matrix. Problem is
%                   that 0 exchange values can be assigned non-zero
%                   radius.
%   * `'max'`       Same as `'min'`, except the positive eigenvalues are
%                   shifted to zero. This works nicely to emphasize
%                   FM type exchange values, has the same problem
%                   as the `'min'` option.
% 
% `'linewidth'`
% : Defines the bond radius if `mode` is set to `line`, one of the folling
%   strings:
%   * `'fix'`       All lines will have a width defined by the `linewidth0`
%                   parameter (default).
%   * `'lin'`       Lines will have a width that is depending 
%                   linearly on the exchange matrix on the bond:
%                          `~Width ~ sum(abs(J))`, 
%                   where the largest line width on
%                   the strongest bond is given by `linewidth0`.
%   * `'pow'`       Same as `'auto'`, but the line width is a
%                   power function of the exchange value:
%                   `W~(sum(abs(J))).^widthpow`.
% 
% `'widthpow'`
% : Defines the power that determines the linewidth if `linewidth`
%   parameter is set to `'pow'`.
% 
% `'linewidth0'`
% : Line width in pt used to draw the bond if `mode` parameter is `'line'`. 
%   Default value is 0.5.
% 
% `'lineStyle'`
% : Determines the line style when `mode` parameter is `'line'`. Possible
%   values are:
%   * `'auto'`      Bonds are plotted as continuous/dashed lines
%                   depending on the label of the corresponding
%                   matrix (dashed line is used if the matrix
%                   label ends with `'-'`, otherwise continuous) (default).
%   * `'--'`        Bonds are plotted as dashed lines.
%   * `'-'`         Bonds are plotted as lines.
% 
% `'zero'`
% : If `true`, bonds with zero exchange matrix will be plotted as
%   well. Default value is `true`.
% 
% `'radius0'`
% : Radius of the cylinder when `mode` parameter is set to `'cylinder'`.
%   Default value is 0.05 \\ang.
% 
% `'radius1'`
% : Radius of the DM vector and the minimum radius of the 
%   ellipsoid, default value is 0.08 \\ang.
% 
% `'radius2'`
% : Constant atom radius, default value is 0.3 \\ang.
% 
% `'radius'`
% : Defines the atom radius (important for arrow bonds, to avoid
%   overlap with the spheres of the atoms), see [swplot.plotatom]:
%   * `'fix'`       Sets the radius of all atoms to the value
%                   given by the `radius2` parameter.
%   * `'auto'`      use radius data from database based on the atom
%                   label multiplied by the `radius2` option value.
% 
% `'ang'`
% : Angle of the arrow head in degree units, default value is 30\\deg.
% 
% `'lHead'`
% : Length of the arrow head, default value is 0.3 \\ang.
% 
% `'scale'`
% : Scaling factor for the length of the DM vector or the size of
%   the ellipsoid relative to the shortest bond length. Default 
%   value is 1/3.
% 
% `'color'`
% : Color of the bonds, one of the following values:
%   * `'auto'`      all bonds gets the color stored in [spinw.matrix],
%   * `'colorname'` all bonds will have the same color defined by the
%                   string, e.g. `'red'`,
%   * `[R G B]`     all bonds will have the same color defined by the given
%                   RGB code.
% 
% `'color2'`
% : Color of the ellipse or DM vector on the bond:
%   * `'auto'`      all object get the color of the bond,
%   * `'colorname'` all object will have the same color defined by the
%                   string, e.g. `'red'`,
%   * `[R G B]`     all objects will have the same color defined by the
%                   given RGB code.
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
% The name of the objects that are created are `'bond'`.
% To find the handles and the corresponding data on these objects, use e.g.
% sObject = swplot.findobj(hFigure,'name','bond')`.
%
% *[DM]: Dzyaloshinskii-Moriya
% *[FM]: FerroMagnet
% *[AFM]: AntiFerromagnet
%

% default values
%fontSize0 = swpref.getpref('fontsize',[]);
nMesh0    = swpref.getpref('nmesh',[]);
nPatch0   = swpref.getpref('npatch',[]);

inpForm.fname  = {'range' 'legend' 'label' 'zero' 'scale' 'radius0' 'mode2' 'linewidth'};
inpForm.defval = {[]      true     true    true   1/3     0.05      []      'fix'      };
inpForm.size   = {[-1 -2] [1 1]    [1 1]   [1 1]  [1 1]   [1 1]     [1 -7]  [1 3]      };
inpForm.soft   = {true    false    false   false  false   false     true    false      };

inpForm.fname  = [inpForm.fname  {'radius' 'mode' 'color' 'nmesh' 'npatch' 'linewidth0' }];
inpForm.defval = [inpForm.defval {'auto'   []     'auto'  nMesh0  nPatch0  0.5          }];
inpForm.size   = [inpForm.size   {[1 -3]   [1 -4] [1 -5]  [1 1]   [1 1]    [1 1]        }];
inpForm.soft   = [inpForm.soft   {false    true  false   false   false    false         }];

inpForm.fname  = [inpForm.fname  {'figure' 'obj' 'unit' 'tooltip' 'radius' 'linestyle' }];
inpForm.defval = [inpForm.defval {[]       []    'lu'   true      'auto'   'auto'      }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1] [1 -6] [1 1]     [1 -8]   [1 -10]     }];
inpForm.soft   = [inpForm.soft   {true     true  false  false     false    false       }];

inpForm.fname  = [inpForm.fname  {'shift' 'replace' 'arrow' 'ang' 'lhead' 'radius2' 'widthpow'}];
inpForm.defval = [inpForm.defval {[0;0;0] true      []      30    0.3     0.3       0.2       }];
inpForm.size   = [inpForm.size   {[3 1]   [1 1]     [1 1]   [1 1] [1 1]   [1 1]     [1 1]     }];
inpForm.soft   = [inpForm.soft   {false   false     true    false false   false     false     }];

inpForm.fname  = [inpForm.fname  {'radius1' 'translate' 'zoom' 'color2' 'copy' 'sign'}];
inpForm.defval = [inpForm.defval {0.08      false        false 'auto'   false  'abs' }];
inpForm.size   = [inpForm.size   {[1 1]     [1 1]        [1 1] [1 -11]  [1 1]  [1 3] }];
inpForm.soft   = [inpForm.soft   {false     false        false false    false  false }];

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
set(hFigure,'Name', 'SpinW: Magnetic bonds');

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
    error('plotbond:WrongInput','The given plotting range is invalid!');
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
        error('plotbond:WrongInput','The given unit string is invalid!');
end

% atom data
mAtom  = obj.matom;

% generate bonds, but don't sort the bonds on DM interactions
[SS, ~] = intmatrix(obj,'plotmode',true,'extend',false,'sortDM',false,'zeroC',param.zero,'nExt',[1 1 1]);

if isempty(SS.all)
    warning('plotbond:EmptyPlot','No bonds to plot!')
    return
end

coupling        = struct;
coupling.dl     = double(SS.all(1:3,:));
coupling.atom1  = SS.all(4,:);
coupling.atom2  = SS.all(5,:);
coupling.matidx = SS.all(15,:);
coupling.idx    = SS.all(16,:);
coupling.cidx   = SS.all(18,:);

% matrix values
mat   = reshape(SS.all(6:14,:),3,3,[]);
% DM interaction
matDM = (mat-permute(mat,[2 1 3]))/2;
% keep symmetric part of the interactions
matSym = (mat+permute(mat,[2 1 3]))/2;
matDM  = permute(cat(2,matDM(2,3,:),matDM(3,1,:),matDM(1,2,:)),[2 3 1]);

% are there non-zero DM vectors
isDM  = any(abs(matDM(:))>5*eps);

matSym0 = obj.matrix.mat(:,:,coupling.matidx);
matSym0 = (matSym0+permute(matSym0,[2 1 3]))/2;
matSym0(2,2,:) = matSym0(2,2,:)-matSym0(1,1,:);
matSym0(3,3,:) = matSym0(3,3,:)-matSym0(1,1,:);
matSym0(1,1,:) = 0;
isSym = any(abs(matSym0(:))>5*eps);

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
    elseif isSym
        param.mode2 = 'sym';
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
mat    = reshape(repmat(permute(mat,[1 2 4 3]),[1 1 nCell 1]),3,3,[]);
matDM  = reshape(repmat(permute(matDM,[1 3 2]),[1 nCell 1]),3,[]);
cIdx   = reshape(repmat(coupling.cidx,[nCell 1]),1,[]);

pos1  = reshape(pos1,3,[]);
pos2  = reshape(pos2,3,[]);

% cut out the bonds that are out of range
switch param.unit
    case 'lu'
        % lower range<=L<= upper range
        pIdx1 = all(bsxfun(@ge,pos1,range(:,1)-10*eps) & bsxfun(@le,pos1,range(:,2)+10*eps),1);
        pIdx2 = all(bsxfun(@ge,pos2,range(:,1)-10*eps) & bsxfun(@le,pos2,range(:,2)+10*eps),1);
        pIdx  = all([pIdx1;pIdx2],1);
    case 'xyz'
        % convert to xyz
        posxyz1 = BV*pos1;
        posxyz2 = BV*pos2;
        pIdx1   = all(bsxfun(@ge,posxyz1,range(:,1)-10*eps) & bsxfun(@le,posxyz1,range(:,2)+10*eps),1);
        pIdx2   = all(bsxfun(@ge,posxyz2,range(:,1)-10*eps) & bsxfun(@le,posxyz2,range(:,2)+10*eps),1);
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
mat    = mat(:,:,pIdx);
matDM  = matDM(:,pIdx);
cIdx   = cIdx(1,pIdx);

% number of bonds to plot
nBond = size(pos1,2);

% bond vector in rlu
dpos = pos2-pos1;
% length of bond in Angstrom
lxyz = sqrt(sum((BV*dpos).^2,1));

% shift positions
pos1 = bsxfun(@plus,pos1,BV\param.shift);
pos2 = bsxfun(@plus,pos2,BV\param.shift);

% color of the bonds
if strcmp(param.color,'auto')
    color = double(obj.matrix.color(:,matidx));
else
    color = repmat(swplot.color(param.color),[1 nBond]);
end

% color of the objects on the bond
if strcmp(param.color2,'auto')
    color2 = color;
else
    color2 = repmat(swplot.color(param.color2),[1 nBond]);
end

% save original matrix values into data
%mat0   = obj.matrix.mat(:,:,matidx);
matDat = [mat;[permute(cIdx,[1 3 2]) zeros(1,2,size(cIdx,2))]];
if numel(matidx)==1
    matDat = mat2cell(matDat,4,3);
else
    matDat = mat2cell(matDat,4,3,ones(1,numel(matidx)));
end

% legend label
lLabel = obj.matrix.label(matidx);

% plot ellipse/arrow on top of bond
switch param.mode2
    case 'none'
        % do nothing, just remove any matrix
        if param.replace
            % find objects to be deleted
            sObj = swplot.findobj(hFigure,'name','bond_mat');
            % delete them!
            swplot.delete([sObj(:).number]);
        end


    case 'antisym'
        % plot arrows
        % remove zero DM vectors
        matDM2 = sum(matDM.^2,1);
        zeroDM = matDM2==0;
        matDM  = matDM(:,~zeroDM);
        posDM  = (pos1(:,~zeroDM)+pos2(:,~zeroDM))/2;
        % maximum value of DM
        maxDM = sqrt(max(matDM2));
        % scale and convert vectors to RLU
        vecDM = BV\((matDM/maxDM)*param.scale*min(lxyz));
        colDM = color2(:,~zeroDM);
        
        swplot.plot('type','arrow','name','bond_mat','position',cat(3,posDM,posDM+vecDM),'text','',...
            'figure',hFigure,'legend',[],'color',colDM,'R',param.radius1,...
            'tooltip',false,'replace',param.replace,'npatch',param.npatch,...
            'data',{},'label',{},'ang',param.ang,...
            'lhead',param.lhead,'translate',param.translate,'zoom',param.zoom);

    case 'sym'
        % plot ellipsoid
        % remove zero ellipsoids
        rmSym = permute(sumn(abs(matSym),[1 2])==0,[1 3 2]);
        % remove Heisenberg interactions
        diagSym = zeros(size(matSym));
        diagSym(1,1,:) = matSym(1,1,:);
        diagSym(2,2,:) = matSym(1,1,:);
        diagSym(3,3,:) = matSym(1,1,:);
        rmSym = rmSym | sum(abs(reshape(bsxfun(@minus,matSym,diagSym),9,[])),1)<10*eps;
        
        matSym = matSym(:,:,~rmSym);
        posSym = (pos1(:,~rmSym)+pos2(:,~rmSym))/2;
        colSym = color2(:,~rmSym);
        % calculating the main radiuses of the ellipsoid.
        [V, Rell] = eigorth(matSym,1e-5);
        % creating positive definite matrix by adding constant to all
        % eigenvalues.
        switch param.sign
            case 'abs'
                % take the absolute value of the exchange
                % problem: this makes AFM and FM couplings equivalent
                Rell = abs(Rell);
            case 'min'
                % shift negative values to zero
                % this makes non-zero from zero values
                Rell = bsxfun(@minus,Rell,min(Rell,[],1).*(min(Rell,[],1)<0));
            case 'max'
                % shift the zero the largest value
                % this makes non-zero from zero values
                Rell = -bsxfun(@minus,Rell,max(Rell,[],1).*(max(Rell,[],1)>0));
            otherwise
                error('plotbond:WrongOption','The given string for the ''sign'' option is invalid!');
        end

        maxR  = sqrt(max(sum(Rell.^2,1)));
        %Rell = ((Rell+param.radius1)/(maxR+param.radius1))*param.scale*min(lxyz);
        Rell = (Rell/maxR)*param.scale*min(lxyz)+param.radius1;
        % V*diag(R) vectorized
        %V = bsxfun(@times,V,permute(Rell,[3 1 2]));
        V = mmat(bsxfun(@times,V,permute(Rell,[3 1 2])),permute(V,[2 1 3]));
        
        swplot.plot('type','ellipsoid','name','bond_mat','position',posSym,'text','',...
            'figure',hFigure,'legend',[],'color',colSym,'T',V,...
            'tooltip',false,'replace',param.replace,'nmesh',param.nmesh,...
            'data',{},'label',{},'nmesh',param.nmesh,'translate',param.translate,...
            'zoom',param.zoom);


    otherwise
        error('plotbond:WrongInput','The given mode2 string is invalid!');
end


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
        pos1 = pos1 + bsxfun(@times,rad1(:)'./lxyz,dpos);
        pos2 = pos2 - bsxfun(@times,rad2(:)'./lxyz,dpos);
    case 'line'
        type0 = 'line';
        % change linestyle based on the matrix label
        % 1 for continuous line
        % 2 for dashed line
        switch param.linestyle
            case 'auto'
                lineStyle0 = (cellfun(@(C)C(end),lLabel)=='-')+1;
            case {'--' '-'}
                lineStyle0 = param.linestyle;
            otherwise
                error('plotbond:WrongInput','The given lineStyle string is illegal!')
        end
    case 'empty'
        type0 = [];
    otherwise
        error('plotbond:WrongInput','The given mode string is illegal!')
end
    
switch param.linewidth
    case 'fix'
        lineWidth = param.linewidth0;
    case 'lin'
        absmat = permute(sumn(abs(mat),[1 2]),[1 3 2]);
        normat = max(absmat);
        if normat<=0
            normat = 1;
        end
        lineWidth = absmat/normat*param.linewidth0;
    case 'pow'
        absmat = permute(sumn(abs(mat),[1 2]),[1 3 2]);
        normat = max(absmat);
        if normat<=0
            normat = 1;
        end
        lineWidth = (absmat/normat).^param.widthpow*param.linewidth0;
    otherwise
        error('plotbond:WrongInput','The given linewidth string is illegal!')
end

% plot bond vectors
if ~isempty(type0)
    swplot.plot('type',type0,'name','bond','position',cat(3,pos1,pos2),'text','',...
        'figure',hFigure,'legend',[],'color',color,'R',param.radius0,...
        'tooltip',false,'replace',param.replace,'lineStyle',lineStyle0,...
        'data',matDat,'label',lLabel,'nmesh',param.nmesh,'ang',param.ang,...
        'lineWidth',lineWidth,'lhead',param.lhead,'translate',param.translate,...
        'zoom',param.zoom,'npatch',param.npatch);
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
    setappdata(hFigure,'legend',lDat);
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

% save range
setappdata(hFigure,'range',struct('range',range,'unit',param.unit));

if param.tooltip
    swplot.tooltip('on',hFigure);
end

end