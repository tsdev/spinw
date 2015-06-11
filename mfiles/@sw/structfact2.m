function sFact = structfact2(obj, hkl, varargin)
% calculates magnetic structure factor using FFT
%
% sFact = STRUCTFACT(obj, k, option1, value1, ...)
%
% Input:
%
% obj       Input sw object, contains positions of the magnetic atoms,
%           size of the magnetic supercell and the vector components of the
%           spins anf g-tensors.
% hkl       Defines the reciprocal lattice vectors where the magnetic
%           intensity is calculated. For commensurate structures these are
%           the possible positions of the magnetic Bragg peaks. For
%           incommensurate helical/conical structures 3 Bragg peaks
%           positions are possible: (k-km,k,k+km) around every reciprocal
%           lattice vector. In this case still the integer positions have
%           to be given and the code calculates the intensities at all
%           three points.
%
% Options:
%
% gtensor   If true, the g-tensor will be included in the static spin
%           correlation function. Including anisotropic g-tensor or
%           different g-tensor for different ions is only possible here.
%           Including a simple isotropic g-tensor is possible afterwards
%           using the sw_instrument() function.
% tol       Tolerance of the incommensurability of the magnetic
%           ordering wavevector. Deviations from integer values of the
%           ordering wavevector smaller than the tolerance are considered
%           to be commensurate. Default value is 1e-4.
%
% Output:
%
% 'sFact' is a structure with the following fields:
% F2        Square of the 3 dimensional magnetic structure factor,
%           dimensions are:
%               [nExt(1)*fExt(1) nExt(2)*fExt(2) nExt(3)*fExt(3)],
%           where nExt is the size of the extended unit cell.
% perp      Square of the perpendicular component of the magnetic structure
%           factor to the Q scattering vector, same size as sFact.int.
% xyz       Cell storing the components of the complex structure factor in
%           the form {X, Y, Z}.
% hkl       Cell storing the wavevectors (h, k, l) in reciprocal lattice
%           units, in the form {h, k, l} where the dimensions of each
%           vector are [1 nExt(ii)*fExt(ii)] where ii is {1, 2, 3} for
%           {h, k, l}.
% obj       Copy of the input obj object.
%
% See also SW_PLOTSF, SW_INTSF, SW.ANNEAL, SW.GENMAGSTR.
%

inpF.fname  = {'gtensor' 'tol'};
inpF.defval = {false     1e-4 };
inpF.size   = {[1 1]     [1 1]};

if iscell(hkl)
    if numel(hkl) ~= 3
        error('sw:structfact:WrongInput','To create a Q grid, give 3 vectors in a cell!');
    end
    % create grid from 3 grid vectors
    [hkl1,hkl2,hkl3] = ndgrid(hkl{:});
    hkl = [hkl1(:);hkl2(:);hkl3(:)];

end
    
% check whether all hkl are integer
if sum(abs(hkl(:)-round(hkl(:))))/numel(hkl) > 1e-10
    error('sw:structfact:WrongInput','All hkl values have to be integers.')
end

param    = sw_readparam(inpF, varargin{:});

matom    = obj.matom;
mAtPos   = matom.r;
nExt     = double(obj.mag_str.N_ext);
nMagAtom = size(mAtPos,2);
km       = obj.mag_str.km;
n        = obj.mag_str.n;

% whether the structure is incommensurate
incomm = any(abs(km-round(km)) > param.tol);

if incomm
    % list of Q-points 
    k = [bsxfun(@minus,hkl,km) hkl bsxfun(@plus,hkl,km)];
end


mMatk = zeros([nExt.*fExt 3]);


% Save magnetic moment directions into the 3D matrices.
mMatX = zeros([nMagAtom nExt]);
mMatY = zeros([nMagAtom nExt]);
mMatZ = zeros([nMagAtom nExt]);

mMatX(:) = obj.mag_str.S(1,:);
mMatY(:) = obj.mag_str.S(2,:);
mMatZ(:) = obj.mag_str.S(3,:);

mMat = reshape(obj.mag_str.S,[3 nExt]);

% Fourier transform + phase in the unit cell.
for ii = 1:nMagAtom
    phiA = exp(-1i*2*pi*sum(hkl.*mAtPos,1));
    
    
    
    mMatk(:,:,:,1) = mMatk(:,:,:,1) + repmat(fftn(shiftdim(mMatX(ii,:,:,:),1)),fExt).*phiA;
    mMatk(:,:,:,2) = mMatk(:,:,:,2) + repmat(fftn(shiftdim(mMatY(ii,:,:,:),1)),fExt).*phiA;
    mMatk(:,:,:,3) = mMatk(:,:,:,3) + repmat(fftn(shiftdim(mMatZ(ii,:,:,:),1)),fExt).*phiA;
    
end


nx  = [0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];
nxn = n'*n;
K1 = 1/2*(eye(3) - nxn - 1i*nx);
K2 = nxn;

% keep the rotation invariant part of Sab
%nx  = [0 -n(3) n(2);n(3) 0 -n(1);-n(2) n(1) 0];
%nxn = n'*n;
m1  = eye(3);

% if the 2*km vector is integer, the magnetic structure is not a true
% helix
%tol = param.tol*2;
%helical =  sum(abs(mod(abs(2*km)+tol,1)-tol).^2) > tol;

if helical
    % integrating out the arbitrary initial phase of the helix
    Sab = 1/2*Sab - 1/2*mmat(mmat(nx,Sab),nx) + 1/2*mmat(mmat(nxn-m1,Sab),nxn) + 1/2*mmat(mmat(nxn,Sab),2*nxn-m1);
end
    

    
% Correction for intensity.
fFact = (nMagAtom*prod(nExt));
mMatk = mMatk/fFact;

% Save the components of the calculated structure factor.
sFact.F2 = sum(abs(mMatk).^2,4);
sFact.xyz = {mMatk(:,:,:,1) mMatk(:,:,:,2) mMatk(:,:,:,3)};

% Save the k-vector values.
sFact.hkl = k';

% Calculate magnetic neutron scattering cross section.
bVect = inv(obj.basisvector);
% Generates Q vectors.
qq{1} = bVect(1,1)*kk{1}+bVect(1,2)*kk{2}+bVect(1,3)*kk{3};
qq{2} = bVect(2,1)*kk{1}+bVect(2,2)*kk{2}+bVect(2,3)*kk{3};
qq{3} = bVect(3,1)*kk{1}+bVect(3,2)*kk{2}+bVect(3,3)*kk{3};
% Normalise Q vectors.
qqNorm = qq{1}.^2+qq{2}.^2+qq{3}.^2;
% Calculates the perpendicular component of F to Q
sFactqq = sFact.xyz{1}.*qq{1}+sFact.xyz{2}.*qq{3}+sFact.xyz{3}.*qq{3};
mMatM = mMatk*0;
mMatM(:,:,:,1) = sFact.xyz{1}-sFactqq.*qq{1}./qqNorm;
mMatM(:,:,:,2) = sFact.xyz{2}-sFactqq.*qq{2}./qqNorm;
mMatM(:,:,:,3) = sFact.xyz{3}-sFactqq.*qq{3}./qqNorm;
sFact.perp = sum(abs(mMatM).^2,4);

sFact.obj = copy(obj);
end