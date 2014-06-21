function [V, rotM] = sw_rot(rotAxis, rotAngle, V)
% rotates vectors around arbitrary axis in 3D
%
% [V, rotM] = SW_ROT(rotAxis, rotAngle, {V})
%
% It rotates vectors in V around rotAxis by rotAngle radian (positive angle
% is the right-hand direction).
%
% Input:
%
% rotAxis   Axis of rotation, dimensions are [1 3].
% rotAngle  Angle of rotation in radian (can be vector with dimensions of
%           [1 nAng]).
% V         Matrix of 3D vectors, dimensions are [3 N], optional.
%
% Output:
%
% V         Rotated vectors, dimensions are [3 N].
% rotM      Rotation matrix, dimensions are [3 3]. If rotAngle is a vector,
%           rotM contains rotation matrices for every angle, it's
%           dimensions are [3 3 nAng].
%
% The rotation matrix defines rotations in a right-handed coordinate
% system, the positive direction is counter-clockwise, when looking from
% where the rotation axis points. To rotate any column vector use the
% following:
%   vp = rotM * v;
%
% To rotate tensors (3x3 matrices) use the following command:
%   Ap = rotM * A * rotM';
%
% See also SW.GENMAGSTR, SW_MIRROR.
%

if nargin==0
    help sw_rot
    return
end

% rotAngle vector along the 3rd dimension
rotAngle = permute(rotAngle(:),[2 3 1]);
% Normalize axis vector.
rotAxis = rotAxis./sqrt(sum(rotAxis.^2));
% 3x3 matrix to calculate Rodrigues' rotation formula.
nx  = [0 -rotAxis(3) rotAxis(2); rotAxis(3) 0 -rotAxis(1); -rotAxis(2) rotAxis(1) 0];
% Rodrigues' rotation formula.
% works for multiple rotation angles as well
% rotM = eye(3)*cos(rotAngle) + sin(rotAngle)*nx + (1-cos(rotAngle))*(rotAxis')*rotAxis;

rotM = bsxfunsym(@times,eye(3),cos(rotAngle)) + bsxfunsym(@times,nx,sin(rotAngle)) + bsxfunsym(@times,(rotAxis')*rotAxis,1-cos(rotAngle));

if nargin > 2
    V = rotM*V;
else
    V = [];
end

end