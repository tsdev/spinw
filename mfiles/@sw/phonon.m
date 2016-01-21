function spectra = phonon(obj, hkl, varargin)
% calculates phonon dispersion based on a spring model
%
% spectra = PHONON(obj, hkl, 'option1', value1 ...)
%
% Input:
%
% obj           Input structure, sw class object.
% hkl           Defines the Q points where the spectra is calculated, in
%               reciprocal lattice units, size is [3 nHkl]. Q can be also
%               defined by several linear scan in reciprocal space. In this
%               case hkl is cell type, where each element of the cell
%               defines a point in Q space. Linear scans are assumed
%               between consecutive points. Also the number of Q points can
%               be specified as a last element, it is 100 by defaults. For
%               example: hkl = {[0 0 0] [1 0 0]  50}, defines a scan along
%               (h,0,0) from 0 to 1 and 50 Q points are calculated along
%               the scan.
%
%               For symbolic calculation at a general reciprocal space
%               point use sym class input. For example to calculate the
%               spectrum along (h,0,0): hkl = [sym('h') 0 0]. To
%               do calculation at a specific point do for example
%               sym([0 1 0]), to calculate the spectrum at (0,1,0).
%
% WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! 
% Works only for Bravais lattice at the moment!!!
% This function is experimental, use it for your own risk!!!
% WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! 
%

% for linear scans create the Q line(s)
if nargin > 1
    hkl = sw_qscan(hkl);
else
    hkl = [];
end


SS = obj.intmatrix('fitmode',2,'conjugate',true);

dR    = SS.all(1:3,:);
%atom1 = SS.all(4,:);
%atom2 = SS.all(5,:);
% average of the diagonals for phonon D
D     = mean(SS.all([6 10 14],:));

% convert dR into xyz coordinate system (Angstrom)
dRxyz = obj.basisvector*dR;
%dRxyz = dR;
dRxyz  = bsxfun(@rdivide,dRxyz,sqrt(sum(dRxyz.^2,1)));

% partial derivative of the distance vector between 2 interacting atoms
% dimensions: [3 3 nBond] --> [alpha beta nBond]
%phiab = bsxfun(@rdivide,bsxfun(@times,permute(dRxyz,[1 3 2]),permute(dRxyz,[3 1 2])),permute(dRl2,[1 3 2]));
phiab = bsxfun(@times,permute(dRxyz,[1 3 2]),permute(dRxyz,[3 1 2]));

%phiab = permute(phiab,[2 1 3]);

% include the spring constant
phiab = bsxfun(@times,permute(D,[1 3 2]),phiab);

% Furier transform exp(-i*2*pi*k*R) factor
% dimensions [1 1 nBond nQ]
dRQ = exp(-1i*2*pi*sum(bsxfun(@times,permute(dR,[4 5 2 3 1]),permute(hkl,[4 5 3 2 1])),5));

% dynamical matrix
% dimensions are [3 3 nQ]
Dab = bsxfun(@plus,-permute(sum(bsxfun(@times,phiab,dRQ),3),[1 2 4 3]),sum(phiab,3));

% solve the eigenvalue problem
[ea,om2] = eigenshuffle(Dab);
om = sqrt(real(om2));

% % X-ray cross section
hklA = (hkl'*2*pi*inv(obj.basisvector))';
int  = permute(sum(bsxfun(@times,ea,permute(hklA,[1 3 2])),1),[2 3 1]).^2;

spectra.Sab   = ea;
spectra.omega = om;
spectra.int   = int;
spectra.hkl   = hkl;
spectra.obj   = copy(obj);
spectra.hklA  = hklA;

end