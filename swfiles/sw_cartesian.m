function [vyOut, vzOut, vxOut] = sw_cartesian(n)
% creates a right handed Cartesian coordinate system
%
% [vy, vz, vx] = SW_CARTESIAN(n)
%
% V = SW_CARTESIAN(n)
%
% It creates an (x,y,z) right handed Cartesian coordinate system with vx,
% vy and vz defines the basis vectors. If only a single output is required,
% V will contain the basis vectors in a matrix: V = [vx vy vz].
%
% Input:
%
% n         Either a 3 element row/column vector or a 3x3 matrix with
%           columns defining 3 vectors.
% Output:
%
% vy,vz,vx  Vectors defining the right handed coordinate system. They are
%           either column of row vectors depending on the shape of the
%           input n.
%


if nargin == 0
    help sw_cartesian
    return
end

% Shape of original vector.
if numel(n) == 3
    nShape = size(n);
    n = n(:);
    
    z = [0; 0;-1];
    y = [0;-1; 0];
    
    if any(cross(n,z))
        vy = cross(n,z);
    else
        vy = cross(n,y);
    end
    vz = cross(n,vy);
    
elseif all(size(n) == [3 3])
    if det(n) == 0
        error('sw_cartesian:WrongInput','The input vectors are not linearly independent!')
    end

    nShape = [3 1];
    vz = cross(n(:,1),n(:,2));
    vy = cross(vz,n(:,1));
    n = n(:,1);
   
else
    error('sw_cartesian:WrongInput','Wrong size of n!')
end


if nargout == 1
    % return a matrix
    vyOut = [n/norm(n) vy/norm(vy) vz/norm(vz)];
else
    % conserve the shape of the input vector.
    vyOut = reshape(vy/norm(vy),nShape);
    vzOut = reshape(vz/norm(vz),nShape);
    vxOut = reshape(n/norm(n),nShape);

end

end