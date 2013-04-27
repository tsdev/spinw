function F2 = structfact(obj, varargin)
% calculates magnetic structure factor using FFT
%
% F2 = STRUCTFACT(swobj, Option1, Value1, ...)
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
% F2.int            Square of the 3 dimensional magnetic structure factor,
%                   dimensions are 
%                   [nExt(1)*fExt(1) nExt(2)*fExt(2) nExt(3)*fExt(3)], 
%                   where nExt is the size of the extended unit cell.
% F2.perp           Square of the perpendicular component of the magnetic
%                   structure factor to the Q scattering vector, same size
%                   as F2.int.
%
% F2.x              X component of the complex structure factor.
% F2.y              Y component of the complex structure factor.
% F2.z              Z component of the complex structure factor.
%
% F2.h              Wavevector h in reciprocal lattice units, dimensions
%                   are [1 nExt(1)*fExt(1)].
% F2.k              Wavevector k in reciprocal lattice units, dimensions 
%                   are [1 nExt(2)*fExt(2)].
% F2.l              Wavevector l in reciprocal lattice units, dimensions 
%                   are [1 nExt(3)*fExt(3)].
% F2.basisvector    Matrix contains the basis vectors in columns.
%
% See also SW_PLOTSF, SW.ANNEAL, SW.GENMAGSTR.
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
F2.int = sum(abs(mMatk).^2,4);
F2.x   = mMatk(:,:,:,1);
F2.y   = mMatk(:,:,:,2);
F2.z   = mMatk(:,:,:,3);

% Save the k-vector values.
F2.h  = 2*pi*k{1};
F2.k  = 2*pi*k{2};
F2.l  = 2*pi*k{3};

% Calculate magnetic neutron scattering cross section.
bVect = inv(obj.basisvector);
% Generates Q vectors.
qq{1} = bVect(1,1)*kk{1}+bVect(1,2)*kk{2}+bVect(1,3)*kk{3};
qq{2} = bVect(2,1)*kk{1}+bVect(2,2)*kk{2}+bVect(2,3)*kk{3};
qq{3} = bVect(3,1)*kk{1}+bVect(3,2)*kk{2}+bVect(3,3)*kk{3};
% Normalise Q vectors.
qqNorm = qq{1}.^2+qq{2}.^2+qq{3}.^2;
% Calculates the perpendicular component of F to Q
F2qq = F2.x.*qq{1}+F2.y.*qq{3}+F2.z.*qq{3};
mMatM = mMatk*0;
mMatM(:,:,:,1) = F2.x-F2qq.*qq{1}./qqNorm;
mMatM(:,:,:,2) = F2.y-F2qq.*qq{2}./qqNorm;
mMatM(:,:,:,3) = F2.z-F2qq.*qq{3}./qqNorm;
F2.perp = sum(abs(mMatM).^2,4);

F2.basisvector = obj.basisvector;
end