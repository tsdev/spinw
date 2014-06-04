function sFact = structfact(obj, varargin)
% calculates magnetic structure factor using FFT
%
% sFact = STRUCTFACT(obj, Option1, Value1, ...)
%
% Input:
%
% obj       Input sw object, contains positions of the magnetic atoms,
%           nExt parameter and the direction of the magnetic spins.
%
% Options:
%
% fExt      Number of reciprocal lattice cells, to calculate form
%           factor, dimensions are [1 3].
%           Default value is [1 1 1].
% S         Optional magnetic structure to use instead of the one stored
%           in obj.
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

inpF.fname  = {'fExt'  'S'           };
inpF.defval = {[1 1 1] obj.mag_str.S };
inpF.size   = {[1 3]   [3 -1]        };

param    = sw_readparam(inpF, varargin{:});

fExt     = param.fExt;
matom    = obj.matom;
mAtPos   = matom.r;
M        = param.S;
nExt     = double(obj.mag_str.N_ext);
nMagAtom = size(mAtPos,2);

mMatk = zeros([nExt.*fExt 3]);

k  = cell(3,1);
kk = cell(3,1);

% Create grid in reciprocal space.
for ii = 1:3
    k{ii} = linspace(0,fExt(ii),nExt(ii)*fExt(ii)+1);
    if numel(k{ii}) == 2
        k{ii} = 0;
    else
        k{ii} = k{ii}(1:end-1);
    end
end

[kk{1}, kk{2}, kk{3}] = ndgrid(k{1},k{2},k{3});

% Save magnetic moment directions into the 3D matrices.
mMatX = zeros([nMagAtom nExt]);
mMatY = zeros([nMagAtom nExt]);
mMatZ = zeros([nMagAtom nExt]);

mMatX(:) = M(1,:);
mMatY(:) = M(2,:);
mMatZ(:) = M(3,:);

% Fourier transform + phase in the unit cell.
for ii = 1:nMagAtom
    
    phiA = exp(-1i*2*pi*(kk{1}*mAtPos(1,ii)+kk{2}*mAtPos(2,ii)+kk{3}*mAtPos(3,ii)));
    mMatk(:,:,:,1) = mMatk(:,:,:,1) + repmat(fftn(shiftdim(mMatX(ii,:,:,:),1)),fExt).*phiA;
    mMatk(:,:,:,2) = mMatk(:,:,:,2) + repmat(fftn(shiftdim(mMatY(ii,:,:,:),1)),fExt).*phiA;
    mMatk(:,:,:,3) = mMatk(:,:,:,3) + repmat(fftn(shiftdim(mMatZ(ii,:,:,:),1)),fExt).*phiA;
    
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