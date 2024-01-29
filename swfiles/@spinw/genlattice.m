function R = genlattice(obj, varargin)
% generates crystal lattice
% 
% ### Syntax
%
% `genlattice(obj,Name,Value)`
%
% `R = genlattice(___)`
%
% ### Description
%
% `genlattice(obj,Name,Value)` generates all necessary parameters to define
% a lattice including space group symmetry and store the result it in the
% [spinw.lattice] field.
%
% `R = genlattice(___)` also returns the rotation matrix that
% rotates the inpub basis vectors to the internal coordinate system.
%
% Alternatively the lattice parameters can be given directly when the
% [spinw] object is created using the `spinw(inpStr)` command, where struct
% contains the fields with initial parameters, e.g.:
% ```
% inpStr.lattice.lat_const = [3 3 4];
% ```
%
% ### Example
%
% ```
% crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'sym','P 6')
% crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'sym',168)
% crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'sym','-y,x-y,z; -x,-y,z','label','R -3 m')
% ```
%
% The three lines are equivalent, both will create hexagonal lattice, with
% $P6$ space group.
%
% ### Input
%
% `obj`
% : [spinw] object.
% 
% ### Options
% 
% `angled`
% : `[\\alpha, \\beta, \\gamma]` angles in \\deg, dimensions are $[1\times 3]$.
% 
% `angle`
% : `[\\alpha, \\beta, \\gamma]` angles in radian, dimensions are $[1\times 3]$.
% 
% `lat_const`
% : `[a, b, c]` lattice parameters in units defined in [spinw.unit] (with \\ang
%   being the default), dimensions are $[1\times 3]$.
% 
% `spgr` or 'sym'
% : Defines the space group. Can have the following values:
%
%   * **space group label** string, name of the space group, can be any
%     label defined in the [symmetry.dat] file.
%   * **space group index** line number in the [symmetry.dat] file.
%   * **space group operators** matrix with dimensions 
%     $[3\times 4\times n_{op}]$.
%   
%   The [symmetry.dat] file stores definition of the 230 space groups in
%   standard settings as it is in the [International Tables of Crystallography](http://it.iucr.org/A/).
%   Additional lines can be added to the [symmetry.dat] file using the
%   [swsym.add] function which later can be used in the `spgr` option.
% 
%   If the `spgr` option is 0, no symmetry will be used. The
%   [spinw.gencoupling] function will determine the equivalent bonds based on
%   bond length.
%   
%   Can also provide spacegroup and label (see below) in a cell e.g.
%   {'-x,y,-z', 'P 2'}
% 
% `label`
% : Optional label for the space group if the generators are given in the
%   `spgr` option.
%
% `bv`
% : Basis vectors given in a matrix with dimensions of $[3\times 3]$, where
%   each column defines a basis vector.
% 
% `origin`
% : Origin for the space group operators, default value is `[0 0 0]`.
% 
% `perm`
% : Permutation of the abc axes of the space group operators.
% 
% `nformula`
% : Gives the number of formula units in the unit cell. It is used
%   to normalize cross section in absolute units. Default value is 0, when
%   cross section is normalized per unit cell.
% 
% ### Output
% 
% `R`
% : Rotation matrix that brings the input basis vector to the SpinW
%   compatible form:
%   ```
%   BVspinw = R*BV
%   ```
% 
% The result of the `spinw.genlattice` function is that `obj.lattice` field
% will be changed based on the input, the lattice parameters are stored
% directly and the optional space group string is converted into space
% group operator matrices.
%
% ### See also
%
% [spinw], [swsym.add], [swsym.operator], [spinw.gencoupling]
%

inpForm.fname  = {'angle'           'lat_const'           'sym'       'label'};
inpForm.defval = {obj.lattice.angle obj.lattice.lat_const []          ''     };
inpForm.size   = {[1 3]             [1 3]                 [-1 -2 -3]  [1 -7] };
inpForm.soft   = {false             false                 true        true   };

inpForm.fname  = [inpForm.fname  {'angled' 'bv'  'spgr'     'origin'           'perm' 'nformula'       }];
inpForm.defval = [inpForm.defval {[0 0 0]  []    []         obj.lattice.origin 'abc'  obj.unit.nformula}];
inpForm.size   = [inpForm.size   {[1 3]    [3 3] [-4 -5 -6] [1 3]              [1 3]  [1 1]            }];
inpForm.soft   = [inpForm.soft   {false    true  true       false              true   true             }];

param = sw_readparam(inpForm, varargin{:});

% input validation
if ~isempty(param.spgr)
    warning('spinw:genlattice:DeprecationWarning',...
            ['spgr parameter name is being deprecated, please use sym',...
            ' instead']);
    if ~isempty(param.sym)
        error('spinw:genlattice:WrongInput', ...
            'Both sym and spgr provided - note sym will be used.');
    else
        param.sym = param.spgr;
    end
end
if any(strcmp('angled', varargin(1:2:end))) && ...
        any(strcmp('angle', varargin(1:2:end)))
	warning('spinw:genlattice:WrongInput', ...
        'Both angle and angled provided - angled will be used.');
end
if isempty(param.sym)
   if norm(param.origin) > 1e-10
      warning('spinw:genlattice:WrongInput', ...
        ['Origin provided without symmetry/spacegroup (both required in ',...
         'same function call) - it will be ignored.']); 
   end
   if ~strcmp(param.perm, 'abc')
     warning('spinw:genlattice:WrongInput', ...
        ['Perm provided without symmetry/spacegroup (both required in ',...
         'same function call) - it will be ignored.']); 
   end
else
    % check valid perm
    invalid_perm_msg = ['Invalid permutation supplied - it must be a ', ...
        'string or numeric array that is a permutation of [1,2,3] or ',...
        'abc respectively.'];
    if ischar(param.perm)
        param.perm = param.perm-'a'+1; % casts to int array e.g. [1,2,3]
    elseif ~isnumeric(param.perm)
        error('spinw:genlattice:WrongInput', invalid_perm_msg);
    end
    % check param.perm is a permutation of [1,2,3]
    if ~any(ismember(perms([1,2,3]), param.perm, 'rows'))
        error('spinw:genlattice:WrongInput', invalid_perm_msg);
    end
    % check valid origin
    if any(param.origin > 1) || any(param.origin < 0)
        error('spinw:genlattice:WrongInput', ...
            'Invalid origin supplied, it must be fractional coordinates.');
    end
end
if iscell(param.sym)
    if numel(param.sym) ~= 2
        error('spinw:genlattice:WrongInput', ...
            ['Cell input for spgr/sym must have two elements {spgr, label}', ...
             ' e.g. {"-x,y,-z", "P 2"}'])
    elseif ~isempty(param.label)
     warning('spinw:genlattice:WrongInput', ...
        ['Label provided in spgr/sym argument and in label argument - ', ...
         'the label will be taken from the label argument']);
    end
end

if ~isempty(param.bv)
    % define basis vector of the new coordinate system
    a = [1 0 0];
    c = [0 0 1];
    
    BV = param.bv;
    % make ab plane perpendicular to z, the starting normal vector:
    n1  = sw_nvect(BV(:,[1 2]),1e-5);
    %n1 = cross(BV(:,1),BV(:,2));
    %n1 = n1/norm(n1);
    % axis of rotation defined by the c x n1
    nRot  = cross(c,n1);
    if norm(nRot) == 0
        nRot = a;
    end
    
    % angle of rotation.
    phi = -atan2(norm(cross(n1,c)),dot(n1,c));
    % rotate the basis vectors
    [BV1, R1] = sw_rot(nRot,phi,BV);
    
    
    % rotate a-axis along x
    phi = atan2(norm(cross(BV1(:,1),a)),dot(BV1(:,1),a));
    % rotate the basis vectors
    if dot(cross(BV1(:,1), BV1(:,2)), BV1(:,3)) > 1e-10
        [BV2, R2] = sw_rot(c, phi, BV1);
    else
        [BV2, R2] = sw_rot(c, -phi, BV1); % ensure BV2 right-handed
    end
    
    
    % check the sign of cross(x,y)
    if c*cross(BV2(:,1),BV2(:,2))<0
        % add extra 180 deg rotation around x
        [BV2, R3] = sw_rot(a,pi,BV2);
        %BV2(:,[2 3]) = -BV2(:,[2 3]);
    else
        R3 = eye(3);
    end
    
    
    % rotation matrix that bring the basis vectors to the right direction
    if nargout > 0
        R = R3*R2*R1;
    end
    
    % the new abc vectors compatible with the SpinW format
    a = BV2(:,1);
    b = BV2(:,2);
    c = BV2(:,3);
    
    obj.lattice.lat_const = [norm(a) norm(b) norm(c)];
    
    angle2 = @(V1,V2)atan2(norm(cross(V1,V2)),dot(V1,V2));
    
    obj.lattice.angle     = [angle2(b,c),angle2(a,c),angle2(a,b)];
    
else
    
    if all(param.angled)
        param.angle = param.angled*pi/180;
    end
    
    obj.lattice.angle     = param.angle;
    obj.lattice.lat_const = param.lat_const;
    
    if nargout > 0
        R = eye(3);
    end
    
end

% copy the apporiate label string
if ~isempty(param.sym) && isempty(param.label)
    if ischar(param.sym)
        param.label = param.sym;
    elseif iscell(param.sym)
        param.label = param.sym{2};
    end
end
    
% generate the symmetry operators
if ~isempty(param.sym)
    if ~iscell(param.sym)
        param.sym = {param.sym};
    end
    [symOp, symInfo] = swsym.operator(param.sym{1});
    % permute the symmetry operators if necessary
    obj.lattice.sym = symOp(param.perm,[param.perm 4],:);
    % assign the origin for space group operators
    obj.lattice.origin = param.origin;
    % set spacegroup label if known from param.sym if the user has not
    % provided a label in the function call
    if  ~isempty(param.label) && ischar(param.label)
        obj.lattice.label = strtrim(param.label);
    elseif isnumeric(param.sym{1}) && numel(param.sym{1})==1
        obj.lattice.label = symInfo.name;
    else
        obj.lattice.label = ''; % can't infer label from matrix of sym ops
    end

else
    if ~isempty(param.label) && ischar(param.label)
        obj.lattice.label = strtrim(param.label);
    end
end   


obj.unit.nformula = int32(param.nformula);

end