function basisVector = basisvector(obj, varargin)
% generates basis vectors of the crystal lattice
%
% basisVec = BASISVECTOR(obj, {norm})
%
% Input:
%
% obj       sw class object.
% norm      If true, the basis vectors are normalized to 1, otherwise the
%           length is equal to the lattice constants. Defult is false.
%           Optional.
%
% Output:
%
% basisVec  Stores the three basis vectors in columns, dimensions are 
%           [3 3].
%
% The 3x3 basisVec matrix can be used also as a coordinate transformation
% matrix from the relative atomic position to positions in the xyz
% coordinate system in Angstrom units.
%
% Example:
%
% To change coordinate system:
%
% relative atomic positions --> xyz
%   r_xyz = basisvector * [ra; rb; rc];
%
% reciprocal lattice units --> Angstrom^-1 (xyz coordinate system)
%   Q_xyz =  [h k l] * 2*pi*inv(basisvector);
%
% See also SW, SW.ABC.
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