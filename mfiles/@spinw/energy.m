function  E = energy(obj, varargin)
% calculates the ground state energy per spin
%
% E = ENERGY(obj, Option1, Value1 ...)
%
% The extended magnetic unit cell, stored in obj, is used for the
% calculation. For non-zero k vector, the interaction energies between
% neighbouring extended unit cells depend on the direction of the moments
% in the two extended unit cells. The angles in further extended unit cells
% are calculated based on the k vector (the k vector is in the units of the
% crystallographic unit cell) and the n vector (normal to the spin rotation
% plane). The moment directions in further extended unit cells are
% calculated by rotating the spins of the extended unit cell in the origin
% by k*R degree around the n vector, where R is the translation vector of
% the origin of the farther extended unit cell. If the extended unit cell
% is equivalent to the crystallographic unit cell, this is equivalent to
% the standard definition of the k vector.
%
% Input:
%
% obj       spinw class object.
%
% Options:
%
% epsilon   The smallest value of incommensurability that is tolerated 
%           without warning. Default is 1e-5.
%
% Output:
%
% E         Energy per moment (anisotropy + exchange + Zeeman energy).
%
%
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
%
% The calculated energy can be wrong for incommensurate structures. For
% example a structure where the spins are rotating in XY plane with an
% incommensurate wavevector of (1/3,0,0). The function only calculates the
% anisotropy energy in the first unit cell, that is for single spin
% Eaniso = Axx*Sxx^2+Ayy*Syy^2. While the anisotropy energy in reality is
% independent of the spin orientation in the XY plane Eaniso=3S*(Axx+Ayy)/2.
% Thus for incommensurate structures one has to be carefull! In the
% triangular case one has to extend the unit cell to nExt = [3 3 1] (in the
% hexagonal setting), in this case the energy will be correct.
%
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
%
% Example:
%
% ...
% cryst.optmagstr('nRun',10)
% E = cryst.energy
%
% After optimising the magnetic structure (by minimizing the ground state 
% energy), the energy per spin is calculated. This can be compared to
% different ground state structures to decide which is the right classical
% ground state of the magnetic model in cryst.
%
% See also SPINW, SPINW.ANNEAL, SPINW.NEWCELL.
%

inpForm.fname  = {'epsilon'};
inpForm.defval = {1e-5     };
inpForm.size   = {[1 1]    };

param = sw_readparam(inpForm,varargin{:});

[SS, SI] = obj.intmatrix;

nExt = double(obj.mag_str.N_ext);
kExt = obj.mag_str.k.*nExt;

% Incommensurate structure in the extended unit cell.
if obj.symbolic
    kTest = sw_sub1(kExt,0.452524243);
else
    kTest = kExt;
end
incomm = any(abs(kTest-round(kTest))>param.epsilon);

M0      = obj.mag_str.S;
nMagExt = size(M0,2);

% Anisotropy energy can be wrong.
if  incomm && (any(any(any((SI.aniso)))) || (any(SI.field)) || (any(any(any(abs(SI.g-SI.g(1)*repmat(eye(3),[1 1 nMagExt]))>param.epsilon)))))
    warning('sw:energy:AnisoFieldIncomm','The calculated energy might be wrong due to incommensurate structure!');
end

if nMagExt>0

    dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
    atom1 = [SS.all(4,:)   1:nMagExt];
    atom2 = [SS.all(5,:)   1:nMagExt];
    JJ    = cat(3,reshape(SS.all(6:14,:),3,3,[]),SI.aniso);
    
    M1 = M0(:,atom1);
    M2 = M0(:,atom2);
    
    % For incommensurate structures in the extended unit cell rotates M2
    % moments.
    if incomm && any(any(dR))
        n = obj.mag_str.n;
        dRIdx = find(any(dR));
        for ii = 1:length(dRIdx)
            M2(:,dRIdx(ii)) = sw_rot(n, kExt*dR(:,dRIdx(ii))*2*pi, M2(:,dRIdx(ii)));
        end
    end
    
    Ml = repmat(permute(M1,[1 3 2]),[1 3 1]);
    Mr = repmat(permute(M2,[3 1 2]),[3 1 1]);
    
    exchE   =  sumsym(sumsym(sumsym(Ml.*JJ.*Mr,3),2),1);
    % TODO
    % correct energy for twins
    ZeemanE = -sumsym(SI.field*permute(mmat(SI.g,permute(M0,[1 3 2])),[1 3 2]),2)*obj.unit.muB;
    
    % energy per spin
    E = (exchE + ZeemanE)/nMagExt;
    
else
    % for undefined magnetic structure, energy returns NaN
    E = NaN;
end

end