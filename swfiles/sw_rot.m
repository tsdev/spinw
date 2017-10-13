function [V, rotM] = sw_rot(rotAxis, rotAngle, V)
% rotates vectorsin 3D
% 
% ### Syntax
% 
% `[~, R] = sw_rot(rotAxis,rotAngle)`
% 
% `[VR, R] = sw_rot(rotAxis,rotAngle,V)`
%
% ### Description
% 
% `[~, R] = sw_rot(rotAxis,rotAngle)` produces the `R` rotation matrix that
% rotates any vector around the given `rotAxis` rotation axis by `rotAngle`
% angle in radian. Positive rotation is the right-hand direction around the
% rotation axis and using the following rotation formula:
% ```
% VR = R*V
% ```
%
% To rotate tensors ($3\times 3$ matrices) use the following formula:
% ```
% Mp = R * M * R';
% ```
%
% `[VR, R] = sw_rot(rotAxis,rotAngle,V)` also rotates the given `V`
% vectors where `VR` are the transformed vectors.
% 
% ### Input Arguments
% 
% `rotAxis`
% : Axis of rotation, stored in a row vector with 3 elements.
% 
% `rotAngle`
% : Angle of rotation in radian, can be also a row vector with $n_{ang}$
%   number of elements.
% 
% `V`
% : Matrix with 3 rows, where each column is a vector in 3D space.
% 
% ### Output Arguments
% 
% `VR`
% : Rotated vectors, stored in a matrix with dimensions of $[3\times N
%   n_{ang}]$.
%
% `R`
% : Rotation matrix with dimensions of $[3\times 3]$ if a single rotation
%   angle is given. If `rotAngle` is a vector, `R` will contain a
%   rotation matrix for each angle, it's dimensions are $[3\times 3\times
%   n_{ang}]$.
%
% ### See Also
% 
% [sw_rotmat] \| [sw_mirror]
%

if nargin==0
    help sw_rot
    return
end

% make row vector
rotAxis = rotAxis(:)';
% rotAngle vector along the 3rd dimension
rotAngle = permute(rotAngle(:),[2 3 1]);
% Normalize axis vector.
rotAxis = rotAxis./sqrt(sum(rotAxis.^2));
% 3x3 matrix to calculate Rodrigues' rotation formula.
nx  = [0 -rotAxis(3) rotAxis(2); rotAxis(3) 0 -rotAxis(1); -rotAxis(2) rotAxis(1) 0];
% Rodrigues' rotation formula.
% works for multiple rotation angles as well
rotM = bsxfunsym(@times,eye(3),cos(rotAngle)) + bsxfunsym(@times,nx,sin(rotAngle)) + bsxfunsym(@times,(rotAxis')*rotAxis,1-cos(rotAngle));

if nargin > 2
    V = mmat(rotM,V);
else
    V = [];
end

end