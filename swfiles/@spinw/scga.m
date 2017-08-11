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
% lambda    If given, the integration is avoided.
% kbase     Basis vectors that span the Brillouin zone if the system is low
%           dimensional. Default value is [] when the dimensionality of the
%           system is determined automatically.
% nQ        Number of Q points where the Brillouin zone is sampled for the
%           integration.
% sublat    List of sublattices.
% isomode   ...
%
% Output:
%
% spectra   Structure with fields:
%   Sab     Spin-spin correlation function stored in a matrix with
%           dimensions of [3,3,D1,D2,...].
%   Dopt    Optimum value of D.
%

% TODO documentation

inpForm.fname  = {'T'    'plot' 'nInt' 'lambda' 'sublat' 'isomode'};
inpForm.defval = {0      true   1e3    []       []       'auto'   };
inpForm.size   = {[1 -1] [1 1]  [1 1]  [1 1]    [1 -2]   [1 -3]   };
inpForm.soft   = {false  false  false  true     true     false    };

param = sw_readparam(inpForm, varargin{:});

if numel(param.T)>1
    param.plot = false;
end

fid  = swpref.getpref('fid',true);

kBT  = param.T*obj.unit.kB;
beta = 1./kBT;
S    = obj.matom.S;

if std(S)~=0 || mean(S)==0
    error('spinw:scga:UnsupportedModel','All magnetic atom has to have the same non-zero spin quantum number!')
end

% number of magnetic atoms in the crystallographic unit cell
if isempty(param.sublat)
    nMag = numel(obj.matom.idx);
else
    nMag = max(param.sublat);
end

fprintf0(fid,'Calculating Self-Consistent Gaussian Approximation (SCGA) (nMag = %d, nQ = %d)...\n',nMag,numel(hkl)/3);

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
    BZ  = sw_qgrid('mat',kbase,'bin',repmat({linspace(0,1,N)},1,D),'fid',0);
    % calculate the Fourier transform of the hamiltonian
    chi = obj.fourier(reshape(BZ,3,[]),'fid',0,'sublat',param.sublat,'isomode',param.isomode);
    % find the eigenvalues over the BZ
    if chi.isiso
        chi.ft = squeeze(chi.ft(1,1,:,:,:));
    else
        chi.ft = reshape(permute(chi.ft,[1 3 2 4 5]),3*nMag,3*nMag,nQBZ);
    end
    
    % number of spin components (Heisenberg=1, otherwise 3)
    nSpinComp = 3-chi.isiso*2;
    
    omega = zeros(nMag*nSpinComp,nQBZ);
    
    for ii = 1:nQBZ
        omega(:,ii) = eig(chi.ft(:,:,ii));
    end
    
    % find the optimum value of lambda
    lambda = zeros(1,numel(beta));
    for ii = 1:numel(beta)
        lambda(ii) = fminsearch(@(lambda)abs(sumn(1./(lambda+beta(ii)*omega),[1 2])/nQBZ/nMag-nSpinComp/3),1-min(beta(ii)*omega(:)));
    end
    
    if param.plot
        % plot the lambda value scan
        nL = 200;
        vL = -min(beta*omega(:))+exp(linspace(-5,log(lambda*1.5+min(beta*omega(:))),nL));
        
        sumA = zeros(1,nL);
        
        for jj = 1:nL
            sumA(jj) = sumn(1./(vL(jj)+beta*omega),[1 2])/nQBZ/nMag;
        end
        
        figure
        loglog(vL,sumA','-')
        line(xlim,nSpinComp/3*[1 1],'color','k')
        hold on
        line(lambda*[1 1],ylim,'color','r')
        xlabel('\lambda(\beta)')
        ylabel('B.Z. sum')
    end
else
    lambda = param.lambda;
end

if numel(beta) == 1
    % calculate spin-spin correlations only if a single temperature is
    % given
    qDim = num2cell(size(hkl));
    
    chi = obj.fourier(reshape(hkl,3,[]),'fid',0,'sublat',param.sublat,'isomode',param.isomode);
    nQ  = numel(hkl)/3;
    
    if chi.isiso
        chi.ft = squeeze(chi.ft(1,1,:,:,:));
    else
        chi.ft = reshape(permute(chi.ft,[1 3 2 4 5]),3*nMag,3*nMag,nQ);
    end
    
    % number of spin components (Heisenberg=1, otherwise 3)
    nSpinComp = 3-chi.isiso*2;
    
    % spin-spin correlation function between any pair of sublattices
    Sabij = bsxfun(@plus,lambda*eye(nMag*nSpinComp),beta*chi.ft);
    
    Sab = zeros(3,3,nQ);
    
    if chi.isiso
        for ii = 1:nQ
            Sab(1,1,ii) = sumn(inv(Sabij(:,:,ii)),[1 2]);
        end
        Sab(3,3,:) = Sab(1,1,:);
        Sab(2,2,:) = Sab(1,1,:);
    else
        for ii = 1:nQ
            Sab(:,:,ii) = squeeze(sumn(reshape(inv(Sabij(:,:,ii)),3,nMag,3,nMag),[2 4]));
        end
    end
    
    % spin-spin sorrelations per magnetic atom
    spec.Sab    = reshape(Sab,3,3,qDim{2:end})/nMag;
    spec.hkl    = hkl;
    spec.hklA   = reshape((reshape(hkl,3,[])'*obj.rl)',3,qDim{2:end});
end

% save lambda only, if multiple temperatures are given
spec.lambda = lambda;
spec.T      = param.T;
spec.isiso  = chi.isiso;

fprintf0(fid,'Calculation finished.\n');

end