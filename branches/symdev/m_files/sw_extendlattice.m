function [mAtom, AAext, SSext] = sw_extendlattice(nExt, mAtom, varargin)
% [mAtom AAext SSext] = SW_EXTENDLATTICE(nExt, mAtom, {AA}, {SS}): produce
% the extended lattice and interactions by tiling the unit cell.
%
% Input:
%
% nExt          Number of unit cell extensions, dimensions are [1 3].
% mAtom         Properties of the magnetic atoms, produced by sw.matom.
% AA            Anisotropy matrices for the magnetic atoms in the unit
%               cell, optional.
% SS            Interactions matrices in the unit cell, optional.
%
% Output:
%
% mAtom         Parameters of the magnetic atoms.
% mAtom.RRext   Positions of magnetic atoms, assuming an extended unit
%               cell, dimensions are [3 nMagExt].
% mAtom.Sext    Spin length of the magnetic atoms, dimensions are
%               [1 nMagExt].
% AAext         Single-ion anisotropy terms in the extended lattice,
%               dimensions are [3 3 nMagExt].
%
% SSext         Interaction matrix in the extended unit cell, struct type.
%               In the struct every field is a matrix. Every column of the
%               matrices describes a single interaction.
% SSext.iso     Isotropic exchange interactions.
% SSext.ani     Anisotropic exchange interations.
% SSext.dm      Dzyaloshinsky-Moriya interaction terms.
% SSext.gen     General 3x3 matrix contains the exchange interaction.
%
% See also SW.INTMATRIX.
%

RR       = mAtom.r;
S        = mAtom.S;
nAtom    = size(RR,2);
nMagAtom = length(S);
nCell    = prod(nExt);
nExt     = nExt';
nExt1    = nExt-1;

rIndex   = 0;
RRext    = zeros(3,nAtom*nCell);
Sext     = zeros(1,nAtom*nCell);

% Extend unit cell.
for kk = 0:nExt1(3)
    for jj = 0:nExt1(2)
        for ii = 0:nExt1(1)
            vIdx = [ii;jj;kk];
            RRext(:,rIndex*nAtom+(1:nAtom)) = (RR + vIdx*ones(1,nAtom))./(nExt*ones(1,nAtom));
            Sext(rIndex*nMagAtom+(1:nMagAtom)) = S;
            rIndex = rIndex + 1;
        end
    end
end

if isfield(mAtom,'idx')
    mAtom.idxext = repmat(mAtom.idx,[1 nCell]);
end

mAtom.RRext = RRext;
mAtom.Sext  = Sext;

switch nargin
    case 2
        AAext = [];
        SSext = struct;
    case 4
        AA = varargin{1};
        SS = varargin{2};
        
        AAext    = zeros(3,3,nMagAtom*nCell);
        fName    = fieldnames(SS);
        SSext    = struct;
        
        for ii = 1:length(fName)
            sName = fName{ii};
            sSize = size(SS.(sName));
            SSext.(sName) = zeros(sSize(1),sSize(2)*nCell);
        end
        
        % Extend anisotropy matrix.
        rIndex = 0;
        for kk = 0:nExt1(3)
            for jj = 0:nExt1(2)
                for ii = 0:nExt1(1)
                    AAext(:,:,rIndex*nMagAtom+(1:nMagAtom)) = AA;
                    rIndex = rIndex + 1;
                end
            end
        end
        
        % Extend coupling matrices.
        for ll = 1:length(fName)
            
            SS0 = double(SS.(fName{ll}));
            N = size(SS0,2);
            SS2 = zeros(size(SS0,1),N*nCell);
            
            if ~isempty(SS0)
                rIndex = 0;
                
                for kk = 0:nExt1(3)
                    for jj = 0:nExt1(2)
                        for ii = 0:nExt1(1)
                            vIdx = [ii;jj;kk];
                            
                            temp = (SS0(1:3,:)+vIdx*ones(1,N))./(nExt*ones(1,N));
                            
                            SS2(1:3,rIndex*N+(1:N))   = floor(temp);
                            SS2(4,rIndex*N+(1:N))     = sum(vIdx.*[1;nExt(1);nExt(1)*nExt(2)],1)*nAtom+SS0(4,:);
                            SS2(5,rIndex*N+(1:N))     = sum((temp-floor(temp)).*(nExt*ones(1,N))*nAtom.*([1;nExt(1);nExt(1)*nExt(2)]*ones(1,N)),1);
                            SS2(5,rIndex*N+(1:N))     = round(SS2(5,rIndex*N+(1:N))+SS0(5,:));
                            SS2(6:end,rIndex*N+(1:N)) = SS0(6:end,:);
                            rIndex = rIndex + 1;
                            
                        end
                    end
                end
            end
            SSext.(fName{ll}) = SS2;
        end
    otherwise
        error('sw:sw_extendunitcell:WrongInput','Wrong number of input parameters!');
end

end