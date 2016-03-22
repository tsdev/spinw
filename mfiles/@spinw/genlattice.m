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
% sym       Space group index, or space group name (string).
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
% crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'sym','P 6')
% crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'sym',168)
%
% The two lines are equivalent, both will create hexagonal lattice, with
% 'P 6' space group.
%
% See also SPINW, SW_ADDSYM, SW_GENSYM, SPINW.GENCOUPLING.
%

inpForm.fname  = {'angle'           'lat_const'           'sym'           'angled' 'bv' };
inpForm.defval = {obj.lattice.angle obj.lattice.lat_const obj.lattice.sym [0 0 0]  []   };
inpForm.size   = {[1 3]             [1 3]                 [1 -1]          [1 3]    [3 3]};
inpForm.soft   = {false             false                 false           false    true };

param = sw_readparam(inpForm, varargin{:});

if ~isempty(param.bv)
    BV = param.bv;
    % make ab plane perpendicular to z
    c = [0 0 1];
    % The starting vector, size (1,3):
    n1 = sw_nvect(BV(:,[1 2]),1e-5);
    % Axis of rotation defined by the spin direction
    nRot  = cross(c,n1);
    % Angle of rotation.
    phi = -atan2(norm(cross(n1,c)),dot(n1,c));
    % rotate the basis vectors
    [BV1, R1] = sw_rot(nRot,phi,BV);
    
    % rotate a-axis along x
    a = [1 0 0];
    phi = -atan2(norm(cross(BV1(:,1),a)),dot(BV1(:,1),a));
    % rotate the basis vectors
    [BV2, R2] = sw_rot(c,-phi,BV1);
    
    % rotation matrix that bring the basis vectors to the right direction
    if nargout > 0
        R = R2*R1;
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

if numel(param.sym) > 1
    % get the number of the symmetry
    [~,~,~,~,param.sym] = sw_gensym(param.sym);
end

obj.lattice.sym       = int32(param.sym);

end