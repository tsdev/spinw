function R = genlattice(obj, varargin)
% generates crystal lattice from given parameters
%
% {R} = GENLATTICE(obj, 'option1', value1, 'option2', value2...)
%
% Input:
%
% obj       spinw class object.
%
% Options:
%
% angled    Alpha, beta, gamma angles in degree, dimensions are [1 3].
% angle     Alpha, beta, gamma angles in radian, dimensions are [1 3].
% lat_const a, b, c lattice parameters, dimensions are [1 3].
% spgr      Space group index, or space group name (string), or space group
%           operators in a matrix with dimensions [3 4 nOp].
% label     Optional label for the space group if the generators are given.
% bv        Basis vectors given in a matrix with dimensions of [3 3].
% origin    Origin for the space group operators. Default is [0 0 0].
% perm      Permutation of the abc axes of the space group operators.
% nformula  Gives the number of formula units in the unit cell. It is used
%           to normalize cross section in absolute units. Default is 0,
%           when cross section is normalized per unit cell.
%
% Output:
%
% R         Rotation matrix that brings the input basis vector to the SpinW
%           compatible form. Optional.
%
% Alternatively the lattice parameters can be given directly when the sw
% object is created using: sw(inpStr), where struct contains the fields
% with initial parameters, e.g.:
%   inpStr.lattice.lat_const = [3 3 4];
%
% The sym option points to the appropriate line in the symmetry.dat file,
% where every line defines a space group by its generators. The first 230
% lines contains all crystallographic space groups with standard setting
% as it is in the International Tables of Crystallography. Additional lines
% can be added to the symmetry.dat file using the sw_addsym() function.
% Every line in the symmetry.dat file can be referenced by either its line
% index or by its label (string).
%
% If the sym option is 0, no symmetry will be used. The sw.gencoupling()
% function will determine the equivalent bonds based on bond length.
%
% Output:
%
% The obj.lattice field will be changed based on the input, the lattice
% constants stored directly and the optional space group string is
% converted to the integer type index.
%
% Example:
%
% ...
% crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'spgr','P 6')
% crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'spgr',168)
%
% The two lines are equivalent, both will create hexagonal lattice, with
% 'P 6' space group.
%
% See also SPINW, SW_ADDSYM, SW_GENSYM, SPINW.GENCOUPLING.
%

inpForm.fname  = {'angle'           'lat_const'           'sym'           'label'          };
inpForm.defval = {obj.lattice.angle obj.lattice.lat_const obj.lattice.sym obj.lattice.label};
inpForm.size   = {[1 3]             [1 3]                 [-1 -2 -3]      [1 -7]           };
inpForm.soft   = {false             false                 false           true             };

inpForm.fname  = [inpForm.fname  {'angled' 'bv'  'spgr'     'origin'           'perm' 'nformula'       }];
inpForm.defval = [inpForm.defval {[0 0 0]  []    []         obj.lattice.origin 'abc'  obj.unit.nformula}];
inpForm.size   = [inpForm.size   {[1 3]    [3 3] [-4 -5 -6] [1 3]              [1 3]  [1 1]            }];
inpForm.soft   = [inpForm.soft   {false    true  true       false              true   true             }];

param = sw_readparam(inpForm, varargin{:});

% new option, but keep the old one as well
if ~isempty(param.spgr)
    param.sym = param.spgr;
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
    [BV2, R2] = sw_rot(c,phi,BV1);
    
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
        %R = R2*R1;
    end
    
    % the new abc vectors compatible with the SpinW format
    a = BV2(:,1);
    b = BV2(:,2);
    c = BV2(:,3);
    
    obj.lattice.lat_const = [norm(a) norm(b) norm(c)];
    obj.lattice.angle     = [sw_angle(b,c),sw_angle(a,c),sw_angle(b,c)];
    
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

% generate the symmetry operators
if ~isempty(param.sym)
    param.sym = sw_gensym(param.sym);
end

% permute the symmetry operators if necessary
if ischar(param.perm)
    param.perm = param.perm-'a'+1;
end
obj.lattice.sym = param.sym(param.perm,[param.perm 4],:);
% assign the origin for space group operators
obj.lattice.origin = param.origin;

if ischar(param.sym)
    param.label = param.sym;
end

obj.lattice.label = param.label;

obj.unit.nformula = int32(param.nformula);

end