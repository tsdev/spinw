function abc = abc(obj,ind)
% returns lattice parameters and angles
%
% latVect = ABC(obj)
%
% Input:
%
% obj       spinw class object
%
% Output:
%
% latVect   Vetor with elements [a, b, c, alpha, beta, gamma],
%           contains the lattice parameters and angles in
%           Angstrom and degree units respectively.
%
% See also SPINW.HORACE.
%

abc = [obj.lattice.lat_const obj.lattice.angle*180/pi];

if nargin>1
    abc = abc(ind);
end

end