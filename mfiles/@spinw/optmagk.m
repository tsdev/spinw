function result = optmagk(obj,varargin)
% determines the magnetic propagation vector
%
% res = OPTMAGK(obj,'option1', value1 ...)
%
% The function determines the optimal propagation vector by calculating the
% Fourier transform of the Hamiltonian as a function of wave vector and
% finding the wave vector that corresponds to the smalles eigenvalue of the
% Hamiltonian. It also returns the normal vector that corresponds to the
% rotating coordinate system. The optimization is achieved via
% Particle-Swarm optimization.
%
% Input:
%
% obj       Input structure, spinw class object.
%
% Options:
%
% Accepts all options of ndbase.pso.
%
% Output:
%
% res       Structure with the following fields:
%               k       Value of the optimal k-vector, with values between 0
%                       and 1/2.
%               n       Normal vector, defines the rotation axis of the
%                       rotating coordinate system.
%               E       The most negative eigenvalue at the given propagation
%                       vector.
%               stat    Full output of the ndbase.pso() optimizer.
%
% See also NDBASE.PSO.
%

%inpForm.fname  = {'fitmode' };
%inpForm.defval = {false     };
%inpForm.size   = {[1 1]     };
%
%param = sw_readparam(inpForm, varargin{:});



% calculate symbolic Fourier transformation if obj is in symbolic mode
if obj.symbolic
    warning('spinw:optmagk:NoSymbolic','The function does not work in symbolic mode!');
    return
end


% generate exchange couplings
[SS, ~, RR] = obj.intmatrix('fitmode',2,'extend',false,'conjugate',true,'zeroC',false);

% list of magnetic atoms in the unit cell
matom = obj.matom;
% number of magnetic atoms in the unit cell
nMagAtom = numel(matom.idx);
% number of bonds
nBond   = size(SS.all,2);

% interacting atom1
%atom1 = SS.all(4,:);
atom1 = mod(SS.all(4,:)-1,numel(matom.S))+1;
% interacting atom1
%atom2 = SS.all(5,:);
atom2 = mod(SS.all(5,:)-1,numel(matom.S))+1;

% exchange energy: J_ij*S_i*S_j
JJ = bsxfun(@times,SS.all(6:14,:),matom.S(atom1).*matom.S(atom2));

% distance vector between interacting atoms in l.u.
dR = SS.all(1:3,:)+RR(:,atom2)-RR(:,atom1);

% number of Q points
nHkl    = 1;

% J^\alpha^\beta component indices in numbers 1:9
idx1 = repmat((1:9)',[1 nBond nHkl]);
% atom1 indices
idx2 = repmat(atom1,[9 1 nHkl]);
% atom2 indices
idx3 = repmat(atom2,[9 1 nHkl]);
% hkl indices
idx4 = repmat(permute(1:nHkl,[1 3 2]),[9 nBond 1]);
% indices all
idxAll = [idx1(:) idx2(:) idx3(:) idx4(:)];

% find the dimensionality of the bonds
% bond vectors of non-zero couplings
dl = SS.all(1:3,:);
% system dimensionality
D = rank(dl);
% k-vector directions
kbase = orth(dl);

kones = ones(1,D);

% optimise the energy using particle swarm
[pOpt0, ~, ~] = ndbase.pso([],@optfun,1/4*kones,'lb',0*kones,'ub',kones,varargin{:},...
    'TolFun',1e-5,'TolX',1e-5,'MaxIter',1e3);

% generate an R-value
optfun2 = @(p)(1e-2+optfun(p)-optfun(pOpt0))^2;

% optimize further using LM
[pOpt, ~, stat] = ndbase.lm2([],optfun2,pOpt0,'TolFun',1e-16,'TolX',1e-6);

kOpt = (kbase*pOpt(:));

% is there any value larger than 1/2
kOpt = mod(kOpt,1);
kOpt(kOpt>1/2) = 1-kOpt(kOpt>1/2);

[Eopt, V] = optfun(pOpt);

% sum up on all atoms
if isreal(V)
    n = [];
else
    %V = sum(reshape(V,3,[]),2);
    V = reshape(V,3,[]);
    % find normal vector
    %n = cross(real(V(1:3)),imag(V(1:3)));
    n = cross(real(V),imag(V));
    n = bsxfun(@rdivide,n,sqrt(sum(n.^2,1)));
    %n = n(:)';
    
    if any(isnan(n))
        warning('spinw:optmagk:NormalVector','The normal vector is undefined, using [001]!')
        n = [];
    end
end

if isempty(n)
    n = repmat([0;0;1],1,nMagAtom);
end

% direction of the imaginary moment
iM = repmat([0;1;0],1,nMagAtom);
rM = cross(iM,n);
nZero = sum(~any(rM,1));
iM(:,~any(rM,1)) = repmat([1;0;0],1,nZero);
rM = cross(iM,n);

% save the optimized values
obj.genmagstr('mode','fourier','k',kOpt','S',1i*iM+rM);

% output results
result.k = kOpt;
result.E = Eopt;
result.n = n;
result.stat = stat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [E, V] = optfun(p)
        %
        
        kvect = (kbase*p(:))';
        
        % exponents, dimension: 1 x nHkl x nBond
        ExpF = permute(exp(1i*2*pi*sum(bsxfun(@times,kvect,permute(dR,[3 1 2])),2)),[2 1 3]);
        
        % J*exp(ikr), dimensions 9 x nBond x nHkl
        Jexp = permute(bsxfun(@times,permute(JJ,[1 3 2]),ExpF),[1 3 2]);
        
        % Exchange energy: A_ij = S_i*S_j*J_ij*exp(ik*2pi*(R_i-R_j))
        % Summed up on all bonds
        % Dimensions: 9 x nMagExt x nMagExt x nHkl
        %                 atom1     atom2
        ft = accumarray(idxAll,Jexp(:)',[9 nMagAtom nMagAtom nHkl]);
        
        if nargout == 1
            % reshape into 3 x 3 x nMagExt x nMagExt x nHkl
            E = min(eig(reshape(permute(reshape(ft,[3 3 nMagAtom nMagAtom nHkl]),[1 3 2 4 5]),3*nMagAtom,3*nMagAtom,[])));
            
        else
            % reshape into 3 x 3 x nMagExt x nMagExt x nHkl
            [V,E] = eig(reshape(permute(reshape(ft,[3 3 nMagAtom nMagAtom nHkl]),[1 3 2 4 5]),3*nMagAtom,3*nMagAtom,[]));
            % find degenerate eigenvalues
            E = diag(E);
            sel = (E-E(1))<10*eps;
            V = sum(V(:,sel),2);
            E = E(1);
        end
    end

end