function rotM = sw_rotmat(rotAxis, rotAngle)
% generates 3D rotation matrix
% 
% ### Syntax
% 
% `R = sw_rotmat(rotAxis,rotAngle)`
%
% ### Description
% 
% `R = sw_rotmat(rotAxis,rotAngle)` produces the `R` rotation matrix that
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
% ### Input Arguments
% 
% `rotAxis`
% : Axis of rotation, stored in a row vector with 3 elements.
% 
% `rotAngle`
% : Angle of rotation in radian, can be also a row vector with $n_{ang}$
%   number of elements.
% 
% ### Output Arguments
% 
% `R`
% : Rotation matrix with dimensions of $[3\times 3]$ if a single rotation
%   angle is given. If `rotAngle` is a vector, `R` will contain a
%   rotation matrix for each angle, it's dimensions are $[3\times 3\times
%   n_{ang}]$.
%
% ### See Also
% 
% [sw_rot] \| [sw_mirror]
%

if nargin==0
    help sw_rotmat
    return
end

[~, rotM] = sw_rot(rotAxis,rotAngle);

end