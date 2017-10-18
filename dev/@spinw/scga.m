function spectra = scga(obj, hkl, varargin)
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
% nInt      Number of Q points where the Brillouin zone is sampled for the
%           integration.
% corr      If true, spinw correlations are calculated at the given
%           momentum vectors. Default is true.
% chi       If true, the magnetic susceptibility is calculated.
% fitmode       Speedup (for fitting mode only), default is false.
% sublat    List of sublattices.
% isomode   ...
% fid       ...
%
% Output:
%
% spectra   Structure with fields:
%   Sab     Spin-spin correlation function stored in a matrix with
%           dimensions of [3,3,1,D1,D2,...].
%   Dopt    Optimum value of D.
%

% TODO documentation

T0   = obj.single_ion.T;

inpForm.fname  = {'T'    'plot' 'nInt' 'lambda' 'sublat' 'isomode' 'fitmode' 'fid'};
inpForm.defval = {T0     true   1e3    []       []       'auto'    false     -1   };
inpForm.size   = {[1 -1] [1 1]  [1 1]  [1 1]    [1 -2]   [1 -3]    [1 1]     [1 1]};
inpForm.soft   = {false  false  false  true     true     false     false     false};

inpForm.fname  = [inpForm.fname  {'formfact' 'formfactfun' 'gtensor' 'corr'        'chi'}];
inpForm.defval = [inpForm.defval {false       @sw_mff      false     ~isempty(hkl) false}];
inpForm.size   = [inpForm.size   {[1 -1]      [1 1]        [1 1]     [1 1]         [1 1]}];
inpForm.soft   = [inpForm.soft   {false       false        false  	 false         false }];

param = sw_readparam(inpForm, varargin{:});

if param.T == 0
    error('spinw:scga:WrongInput','Invalid temperature, the SCGA solver needs non-zero temperature!')
end

if ~param.fitmode
    % save the time of the beginning of the calculation
    spectra.datestart = datestr(now);
end

if param.fid == -1
    fid  = swpref.getpref('fid',true);
else
    fid = param.fid;
end


if numel(param.T)>1
    param.plot = false;
end

kBT  = param.T*obj.unit.kB;
beta = 1./kBT;
nT   = numel(beta);
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
    SS = obj.intmatrix(struct('fitmode',true,'extend',false,'conjugate',true,'zeroC',false),'noCheck');
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
    chiFT = obj.fourier(reshape(BZ,3,[]),struct('fid',0,'sublat',param.sublat,'isomode',param.isomode),'noCheck');
    % find the eigenvalues over the BZ
    if chiFT.isiso
        chiFT.ft = permute(chiFT.ft(1,1,:,:,:),[3 4 5 1 2]);
    else
        chiFT.ft = reshape(permute(chiFT.ft,[1 3 2 4 5]),3*nMag,3*nMag,nQBZ);
    end
    
    % TODO
    %chiFT.ft = bsxfun(@plus,chiFT.ft,eye(size(chiFT.ft,1)));
    % number of spin components (Heisenberg=1, otherwise 3)
    nSpinComp = 3-chiFT.isiso*2;
    
    omega = zeros(nMag*nSpinComp,nQBZ);
    
    for ii = 1:nQBZ
        omega(:,ii) = eig(chiFT.ft(:,:,ii));
    end
    
    % find the optimum value of lambda
    lambda = zeros(1,numel(beta));
    for ii = 1:numel(beta)
        lambda(ii) = fminsearch(@(lambda)abs(sumn(1./(lambda+beta(ii)*omega),[1 2])/nQBZ/nMag-nSpinComp/3),1-min(beta(ii)*omega(:)));
    end
    
    if param.chi
        chi =  bsxfun(@plus,bsxfun(@times,reshape(repmat(lambda,[nQBZ 1]),[1 1 nQBZ*nT]),eye(nMag*nSpinComp)),...
            bsxfun(@times,reshape(repmat(beta,  [nQBZ 1]),[1 1 nQBZ*nT]),repmat(chiFT.ft,[1 1 nT])));
        if chiFT.isiso
            chi = sum(reshape(3*sumn(invfast(chi).^2,[1 2]),[nQBZ nT]))/nMag.*beta;
        else
            chi = sum(reshape(sumn(invfast(chi).^2,[1 2]),[nQBZ nT]))/nMag.*beta;
        end
        spectra.chi = chi;
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

if param.corr
    % calculate spin-spin correlations
    qDim        = size(hkl);
    spectra.hkl = hkl;
    hkl         = reshape(hkl,3,[]);
    chiFT         = obj.fourier(hkl,struct('fid',0,'sublat',param.sublat,'isomode',param.isomode),'noCheck');
    nQ          = numel(hkl)/3;
    
    if chiFT.isiso
        chiFT.ft = permute(chiFT.ft(1,1,:,:,:),[3 4 5 1 2]);
    else
        chiFT.ft = reshape(permute(chiFT.ft,[1 3 2 4 5]),3*nMag,3*nMag,nQ);
    end
    
    % number of spin components (Heisenberg=1, otherwise 3)
    nSpinComp = 3-chiFT.isiso*2;
    
    % spin-spin correlation function between any pair of sublattices
    % Sabij: [nDegFree nDegFree nHkl*nT]
    %
    % Sabij = bsxfun(@plus,permute(lambda*eye(nSub),[3 4 1 2]),beta*chi.ft);
    Sabij =  bsxfun(@plus,bsxfun(@times,reshape(repmat(lambda,[nQ 1]),[1 1 nQ*nT]),eye(nMag*nSpinComp)),...
                          bsxfun(@times,reshape(repmat(beta,  [nQ 1]),[1 1 nQ*nT]),repmat(chiFT.ft,[1 1 nT])));
    
    if chiFT.isiso
        Sab = permute(invfast(Sabij),[4 5 1:3]);
        Sab(3,3,:,:,:) = Sab(1,1,:,:,:);
        Sab(2,2,:,:,:) = Sab(1,1,:,:,:);
    else
        Sab = permute(reshape(invfast(Sabij),3,nMag,3,nMag,[]),[1 3 2 4 5]);
    end
    
    if param.formfact
        spectra.formfact = true;
        % Angstrom^-1 units for Q
        hklA = (hkl'*obj.rl)';
        % store form factor per Q point for each atom in the magnetic supercell
        FF = param.formfactfun(permute(obj.unit_cell.ff(1,:,obj.matom.idx),[3 2 1]),hklA);
        % include the form factor in the z^alpha, z^beta matrices
        FF  = bsxfun(@times,permute(FF,[3 4 5 1 2]),permute(FF,[3 4 1 5 2]));
        Sab = bsxfun(@times,repmat(FF,[1 1 1 1 nT]),Sab);
        spectra.hklA   = reshape(hklA,[3,qDim(2:end)]);
    else
        spectra.formfact = false;
        spectra.hklA     = reshape((hkl'*obj.rl)',[3,qDim(2:end)]);
    end

    % spin-spin sorrelations per magnetic atom
    spectra.Sab    = reshape(permute(sumn(Sab,[3 4]),[1 2 5 3 4]),[3,3,1,qDim(2:end),nT])/nMag;
    % save the new temperature
    if nT == 1
        obj.single_ion.T = param.T;
    end
end

% save lambda only, if multiple temperatures are given
spectra.lambda   = lambda;
spectra.T        = param.T;
spectra.isiso    = chiFT.isiso;
spectra.formfact = param.formfact;
spectra.gtensor  = param.gtensor;

if ~param.fitmode
    spectra.obj     = copy(obj);
    spectra.dateend = datestr(now);
end

fprintf0(fid,'Calculation finished.\n');

end