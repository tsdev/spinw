function optRes = tisza(obj,varargin)
% find magnetic structure using Luttinger-Tisza method
%
% optRes = TISZA(obj,'option1', value1 ...)
%
% The function determines the magnetic structure using the Luttinger-Tisza
% method.
%
% Options:
%
% n         Normal vector of the spiral plane, default is [0 0 1].
% hkl       Points in momentum space for minimizing the energy.
% sub       Row vector of indices, that assign each magnetic moment to a
%           sublattice. The length of the vector is the number of magnetic
%           atoms in the magnetic supercell. Default value is 1:nMagExt,
%           where each magnetic moment belongs to a different sublattice.
%
% See also SPINW.OPTMAGSTR, SPINW.OPTMAGSTEEP.
%

% magnetic atoms in the unit cell
matom = obj.matom;

% number of magnetic atoms in the magnetic supercell
nMagExt = prod(obj.mag_str.N_ext)*numel(matom.idx);

inpForm.fname  = {'n'     'sub'     };
inpForm.defval = {[0 0 1] 1:nMagExt };
inpForm.size   = {[1 3]   [1 -1]    };

param = sw_readparam(inpForm, varargin{:});

% Fourier transform of the magnetic interactions
if obj.symbolic
    % dimensions: 3 x 3 x nMagExt x nMagExt
    ftRes = fourier(obj);
    ft = ftRes.ft;
else
	% dimensions: 3 x 3 x nMagExt x nMagExt x nHkl
    % ...
    return
end

% sublattice vector
sub = param.sub;
% number of moments in each sublattice
subC = histcounts(sub);
% sublattice index that needs summation
subSumIdx = find(subC>1);
% number of sublattices
nSub = max(sub);

if ~all(subC) || sum(subC)~=nMagExt
    error('sw:tisza:WrongInput','Wrong list of sublattices!')
end

% reduce the nMagExt x nMagext matrix to nSub * nSub
for ii = 1:numel(subSumIdx)
    % sum up rows of the Fourier transform that belongs to the same
    % sublattice
    sumIdx = find(sub ==subSumIdx(ii));
    ft(:,:,subSumIdx(ii),:,:) = sumsym(ft(:,:,sumIdx,:,:),3);
    % remove the summed up rows
    ft(:,:,sumIdx(2:end),:,:) = [];
    % sum up columns of the Fourier transform that belongs to the same
    % sublattice
    ft(:,:,:,subSumIdx(ii),:) = sumsym(ft(:,:,:,sumIdx,:),4);
    % remove the summed up columns
    ft(:,:,:,sumIdx(2:end),:) = [];
    % also remove the deleted elements from the sublattice vector
    sub(sumIdx(2:end)) = [];
end

optRes.ft  = ft;
%optRes.hkl = 
end