function basisVector = basisvector(obj)
% generates basis vectors
%
% basisVector = BASISVECTOR(obj)
%
% basisVector   Stores the three basis vectors in columns, dimensions are
%               [3 3].
%

alpha = obj.lattice.angle(1);
beta  = obj.lattice.angle(2);
gamma = obj.lattice.angle(3);

v1 = [1          0          0];
v2 = [cos(gamma) sin(gamma) 0];

v3(1) = cos(beta);
v3(2) = sin(beta)*(cos(alpha)-cos(beta)*cos(gamma))/(sin(beta)*sin(gamma));
v3(3) = sqrt(sin(beta)^2-v3(2)^2);

basisVector = [v1' v2' v3']*diag(obj.lattice.lat_const);

end