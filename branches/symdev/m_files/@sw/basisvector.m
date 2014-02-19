function basisVector = basisvector(obj, varargin)
% generates basis vectors
%
% basisVector = BASISVECTOR(obj, {norm})
%
% basisVector   Stores the three basis vectors in columns, dimensions are
%               [3 3].
%
% basisVector can be used also as a coordinate transformation matrix. If
% norm is true, the basis vectors will be normalized. Default is false.
%
% To change coordinate system:
%
% relative atomic positions --> xyz
%   r_xyz = basisvector * [ra; rb; rc];
%
% reciprocal lattice units --> Angstrom^-1 (xyz coordinate system)
%   Q_xyz =  [h k l] * 2*pi*inv(basisvector);
%

if nargin == 1
    norm = false;
else
    norm = varargin{1};
end

alpha = obj.lattice.angle(1);
beta  = obj.lattice.angle(2);
gamma = obj.lattice.angle(3);

v1 = [1          0          0];
v2 = [cos(gamma) sin(gamma) 0];

v3(1) = cos(beta);
v3(2) = sin(beta)*(cos(alpha)-cos(beta)*cos(gamma))/(sin(beta)*sin(gamma));
v3(3) = sqrt(sin(beta)^2-v3(2)^2);

basisVector = [v1' v2' v3'];

% if not normalized, use the lattice constants
if ~norm
    basisVector = basisVector*diag(obj.lattice.lat_const);
end

end