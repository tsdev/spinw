function abc = abc(obj)
% returns lattice parameters and angles
%
% latVect = ABC(obj)
%
% Input:
%
% obj       sw class object
%
% Output:
%
% latVect   Vetor with elements [a, b, c, alpha, beta, gamma],
%           contains the lattice parameters and angles in
%           Angstrom and degree units respectively.
%
% See also SW.HORACE.
%

abc = [obj.lattice.lat_const obj.lattice.angle*180/pi];

end