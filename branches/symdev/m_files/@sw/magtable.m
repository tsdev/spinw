function moments = magtable(obj)
% creates tabulated list of all magnetic moments stored in sw object
%
% moments = MAGTABLE(obj)
%
% The function lists the moment directions in the magnetic supercell, whose
% size is defined by the obj.mag_str.N_ext field. The positions of the
% magnetic atoms are in lattice units.
%
% Output:
%
% Struct tye data that contains the following fields:
%   M       Matrix, where every columndefines a magnetic moment, dimensions
%           are [3 nMagExt].
%   R       Matrix, where evry column defines the position of the magnetic
%           atom in lattice units.
%
% See also SW.GENMAGSTR.
%

moments.M = obj.mag_str.S;

mAtom = obj.matom;
nExt  = double(obj.mag_str.N_ext);

% Create mAtom.Sext matrix.
mAtom    = sw_extendlattice(nExt, mAtom);

if isempty(moments.M)
    moments.R = zeros(3,0);
else
    % Positions in l.u.
    moments.R = bsxfun(@times,mAtom.RRext,nExt');
end

end