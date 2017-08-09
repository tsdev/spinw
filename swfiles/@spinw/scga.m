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

inpForm.fname  = {'T'   'plot' 'Dlim'  'kbase' 'nQ'  'Dopt'};
inpForm.defval = {0     true   [0 0.5] []      1e3   []    };
inpForm.size   = {[1 1] [1 1]  [1 2]   [3 -4]  [1 1] [1 1] };
inpForm.soft   = {false false  false   true    false true  };

param = sw_readparam(inpForm, varargin{:});

kBT  = sw_converter(param.T,obj.unit.label{4},obj.unit.label{2});
beta = 1/kBT;
S    = obj.unit_cell.S(obj.unit_cell.S>0);

if std(S)~=0 || mean(S)==0
    error('spinw:scga:UnsupportedModel','All magnetic atom has to have the same non-zero spin quantum number!')
else
    S = S(1);
end

if isempty(param.Dopt)
    if isempty(param.kbase)
        % find the dimensionality of the bonds
        % bond vectors of non-zero couplings
        % generate exchange couplings
        SS = obj.intmatrix('fitmode',2,'extend',false,'conjugate',true,'zeroC',false);
        dl = SS.all(1:3,:);
        % system dimensionality
        D = rank(dl);
        % k-vector directions
        kbase = orth(dl);
    else
        kbase = param.kbase;
        D = size(kbase,2);
    end
    
    % q-points
    N    = round(param.nQ^(1/D));
    BZ  = sw_qgrid('mat',kbase,'bin',repmat({linspace(0,1,N)},1,D));
    
    chi0 = obj.fourier(reshape(BZ,3,[]));
    FT   = chi0.ft;
    
    % optimised solution
    Dopt = sw_fminsearchbnd(@(D)saddle_eq(D,FT,beta,S),mean(param.Dlim),param.Dlim(1),param.Dlim(2));
    
    if param.plot
        nD = 200;
        Dv = linspace(param.Dlim(1),param.Dlim(2),nD);
        
        sumA = zeros(1,nD);
        
        nMat = size(FT,5);
        
        for jj = 1:nD
            sumA(jj) = sumn(ndbase.inv3(bsxfun(@plus,FT,Dv(jj)*eye(3)),'diag'),1:4)/nMat/beta;
        end
        
        figure
        plot(Dv,sumA','o-')
        line(xlim,S^2*[1 1],'color','k')
        hold on
        line(Dopt*[1 1],ylim,'color','r')
        xlabel('\Delta(\beta)')
        ylabel('B.Z. sum')
    end
else
    Dopt = param.Dopt;
end

% calculate spin-spin correlations
qDim = num2cell(size(hkl));

chi = obj.fourier(reshape(hkl,3,[]));

A0   = beta*bsxfun(@plus,chi.ft,Dopt*eye(3));

spec.Sab  = reshape(permute(sumn(real(ndbase.inv3(A0)),[3 4]),[1 2 5 3 4]),3,3,qDim{2:end});
spec.Dopt = Dopt;
spec.hkl  = hkl;
spec.hklA = reshape((reshape(hkl,3,[])'*obj.rl)',3,qDim{2:end});

end

function R = saddle_eq(D,FT,beta,S)
% saddle_eq(D,FT,beta)
% to solve saddle point equation

nMat = size(FT,5);

sumA = sumn(ndbase.inv3(bsxfun(@plus,FT,D*eye(3)),'diag'),1:5)/nMat/beta;

R = abs(real(sumA-S^2));

end