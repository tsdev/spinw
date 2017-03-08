function phi = sw_angle(V1,V2)
% calculates the angle between 2 vectors
%
% phi = SW_ANGLE(V1,V2)
%
% Input:
%
% V1    Input vector with 3 elements.
% V2    Input vector with 3 elements.
%
% Output:
%
% phi   Angle in radian.
%

if nargin == 0
    help sw_angle;
    return
end

phi = atan2(norm(cross(V1,V2)),dot(V1,V2));

end