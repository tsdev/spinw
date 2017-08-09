function spec = scga(obj, hkl, varargin)
% applies the self consistent Gaussian approximation at finite temperature
%
% spectra = SCGA(obj, hkl, 'option1', value1 ...)
%
% Input:
%
% obj       SpinW object.
% hkl       Defines the Q points where the correlations are calculated. It
%           is a matrix with dimensions [3,D1,D2,...]. Where the first
%           dimesion corresponds to the [h,k,l] index of the Q-point.
%
% Options:
%
% T         Temperature of the calculation in units given by obj.unit.
% plot      If true, the fitting of the integration constant is plotted.
% Dlim      Limits of the integration constant.
% Dopt      If given, the integration is avoided.
% kbase     Basis vectors that span the Brillouin zone if the system is low
%           dimensional. Default value is [] when the dimensionality of the
%           system is determined automatically.
% nQ        Number of Q points where the Brillouin zone is sampled for the
%           integration.
%
% Output:
%
% spectra   Structure with fields:
%   Sab     Spin-spin correlation function stored in a matrix with
%           dimensions of [3,3,D1,D2,...].
%   Dopt    Optimum value of D.
%

inpForm.fname  = {'T'   'plot' 'nInt' 'lambda' 'subLat'};
inpForm.defval = {0     true   1e3    []       []      };
inpForm.size   = {[1 1] [1 1]  [1 1]  [1 1]    [1 -1]  };
inpForm.soft   = {false false  false  true     true    };

param = sw_readparam(inpForm, varargin{:});

kBT  = param.T*obj.unit.kB;
beta = 1/kBT;
S    = obj.matom.S;

if std(S)~=0 || mean(S)==0
    error('spinw:scga:UnsupportedModel','All magnetic atom has to have the same non-zero spin quantum number!')
end


% number of magnetic atoms in the crystallographic unit cell
nMag = numel(obj.matom.idx);

% sublattices
subLat = param.subLat;

% reduce the lattice into sublattices if necessary
if ~isempty(subLat)
    % number of sublattices
    nSub = max(subLat);
    idx1 = repmat(subLat(:),1,nMag);
    idx2 = repmat(subLat(:)',nMag,1);
    subs = [idx1(:) idx2(:)];
else
    nSub = nMag;
end

if isempty(param.lambda)
    % find the dimensionality of the bonds
    % bond vectors of non-zero couplings
    % generate exchange couplings
    SS = obj.intmatrix('fitmode',true,'extend',false,'conjugate',true,'zeroC',false);
    % calculate the basis vectors
    L = sw_bonddim(SS.all(1:5,:));
    % unite all basis vectors to get the dimensionality of the full system
    dl    = [L(:).base];
    kbase = orth(dl);
    % system dimensionality
    D     = size(kbase,2);
    
    % q-points
    N    = round(param.nInt^(1/D));
    nQBZ = N^D;
    BZ  = sw_qgrid('mat',kbase,'bin',repmat({linspace(0,1,N)},1,D));
    
    chi0 = obj.fourier(reshape(BZ,3,[]));
    % include the spin value into the Fourier transform of the Js
    % thus convert the model into interacting S=1 spins
    FT = bsxfun(@times,chi0.ft,permute(bsxfun(@times,S',S),[3 4 1 2]));

    % reduce the lattice into sublattices
    if ~isempty(subLat)
        FT  = reshape(FT,3,3,[],nQBZ);
        FT2 = zeros(3,3,nSub,nSub,nQBZ);
        
        for ii = 1:nSub
            for jj = 1:nSub
                FT2(:,:,ii,jj,:) = permute(sum(FT(:,:,ismember(subs,[ii jj],'rows'),:),3),[1 2 3 5 4])*nSub/nMag;
                if ii == jj
                    FT2(:,:,ii,ii,:) = FT2(:,:,ii,ii,:) + 1;
                end
            end
        end
        FT = FT2;
    else
        for ii = 1:nSub
            FT(:,:,ii,ii,:) = FT(:,:,ii,ii,:) + 1;
        end
    end
    % find the eigenvalues over the BZ
    [~,omega] = eigorth(squeeze(FT(1,1,:,:,:)));

    % find the optimum value of lambda
    lambda = fminsearch(@(lambda)abs(sumn(1./(lambda+beta*omega),[1 2])/nQBZ/4-1/3),3)
    
    
    if param.plot
        % plot the lambda value scan
        nL = 200;
        vL = linspace(0,3,nL);
        
        sumA = zeros(1,nL);
        
        for jj = 1:nL
            sumA(jj) = sumn(1./(vL(jj)+beta*omega),[1 2])/nQBZ*nSub/nMag;
        end
        
        figure
        plot(vL,sumA','-')
        line(xlim,[1/3 1/3],'color','k')
        hold on
        line(lambda*[1 1],ylim,'color','r')
        xlabel('\lambda(\beta)')
        ylabel('B.Z. sum')
    end
else
    lambda = param.lambda;
end

% calculate spin-spin correlations
qDim = num2cell(size(hkl));

chi = obj.fourier(reshape(hkl,3,[]));
nQ  = numel(hkl)/3;

FT = chi.ft;
% include the spin value into the Fourier transform of the Js
% thus convert the model into interacting S=1 spins
FT = bsxfun(@times,FT,permute(bsxfun(@times,S',S),[3 4 1 2]));

% reduce the lattice into sublattices
if ~isempty(subLat)
    FT  = reshape(FT,3,3,[],nQ);
    FT2 = zeros(3,3,nSub,nSub,nQ);
    
    for ii = 1:nSub
        for jj = 1:nSub
            FT2(:,:,ii,jj,:) = permute(sum(FT(:,:,ismember(subs,[ii jj],'rows'),:),3),[1 2 3 5 4])*nSub/nMag;
            if ii == jj
                FT2(:,:,ii,ii,:) = FT2(:,:,ii,ii,:) + 1;
            end
        end
    end
    FT = FT2;
else
    for ii = 1:nSub
        FT(:,:,ii,ii,:) = FT(:,:,ii,ii,:) + 1;
    end
end

% spin-spin correlation function between any pair of sublattices
Sabij = bsxfun(@plus,permute(lambda*eye(nSub),[3 4 1 2]),beta*FT);

Sab = zeros(1,nQ);
for ii = 1:nQ
    Sab(ii) = sumn(inv(squeeze(Sabij(1,1,:,:,ii))),[1 2]);
end

% spin-spin sorrelations per magnetic atom
spec.Sab    = reshape(Sab,qDim{2:end})/nSub;
spec.lambda = lambda;
spec.hkl    = hkl;
spec.hklA   = reshape((reshape(hkl,3,[])'*obj.rl)',3,qDim{2:end});

end