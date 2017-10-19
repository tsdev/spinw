function rotM = sw_rotmatd(rotAxis, rotAngle)
% generates 3D rotation matrix
% 
% ### Syntax
% 
% `R = sw_rotmatd(rotAxis,rotAngle)`
%
% ### Description
% 
% `R = sw_rotmatd(rotAxis,rotAngle)` produces the `R` rotation matrix, for
% identically to [sw_rotmat], except that here the unit of `rotAngle` is
% \\deg.
%
% ### See Also
%
% [sw_rotmat] \| [sw_rot]
%

if nargin==0
    help sw_rotmatd
    return
end

[~, rotM] = sw_rot(rotAxis,rotAngle*pi/180);

end