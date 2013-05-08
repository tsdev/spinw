function gencoupling(obj, varargin)
% generates the COUPLING property of sw object
%
% GENCOUPLING(obj, 'option1', value, ...)
%
% It calculates the symmetry equivalent couplings between magnetic atoms
% and sorts them according to the atom-atom distance. If option 'sym' is
% false, the crystal symmetry is not considered.
%
% Options:
%
% sym           If true, simmetry equivalent pairs of atoms are calculated.
%               If false, equivalent pairs will be identified based on
%               distance. If the crystal has space group symmetry > 1, the
%               default is true, otherwise false.
% nUnitCell     Edge length of the parallelepiped withing the
%               algorithm searches for neares neighbours in lattice
%               units. Default is 3.
% maxDistance   Maximum inter-ion distance that will be stored in the
%               obj.coupling property in units of Angstrom. Default
%               is 6.
% tol           Tolerance of distance, within two couplings are regarded
%               equivalent, default is 1e-5.
%
% See also SW.
%

isSym = obj.lattice.sym > 1;

inpForm.fname  = {'sym' 'nUnitCell' 'maxDistance' 'tol' };
inpForm.defval = {isSym 3           6             1e-5  };
inpForm.size   = {[1 1] [1 1]       [1 1]         [1 1] };

param = sw_readparam(inpForm, varargin{:});
tol   = param.tol;

% Number of unit cells along any direction to find nearest neighbours.
nUnitCell = param.nUnitCell;

% Generate all atomic positions in the unit cell.
mAtom      = obj.matom;

% Number of magnetic atoms in the unit cell.
nMagAtom = size(mAtom.r,2);

if nMagAtom > 0
    
    % Define number of magnetic atoms in the half 'cube' around the center.
    nHalfCube = (4*nUnitCell^3+6*nUnitCell^2+3*nUnitCell)*nMagAtom;
    
    % Unit cell translation
    atTr    = zeros(3,nHalfCube);
    % Atomic position
    atPos   = zeros(3,nHalfCube);
    % Atom index
    atIndex = zeros(1,nHalfCube);
    
    
    % Get all the magnetic atomic positions of in the half 'cube'.
    index = 1;
    for ii = 1:nUnitCell
        for jj = 1:(2*nUnitCell+1)
            for kk = 1:(2*nUnitCell+1)
                for ll = 1:nMagAtom
                    atTr(:,index)  = [ii; (jj-nUnitCell-1); (kk-nUnitCell-1)];
                    atPos(:,index) = mAtom.r(:,ll) + atTr(:,index);
                    atIndex(index) = ll;
                    index = index+1;
                end
            end
        end
    end
    for jj = 0:nUnitCell
        for kk = ((jj==0)*(nUnitCell+1)+1):(2*nUnitCell+1)
            for ll = 1:nMagAtom
                atTr(:,index)  = [0; jj; (kk-nUnitCell-1)];
                atPos(:,index) = mAtom.r(:,ll) + atTr(:,index);
                atIndex(index) = ll;
                index = index+1;
            end
        end
    end
    
    % Number of elements in the neighbour list.
    nDist = size(atIndex,2)*nMagAtom + nMagAtom*(nMagAtom-1)/2;
    % matrix stores [dlx dly dlz matom1 matom2 idx atom1 atom2 cx cy cz dist dx dy dz]
    sortM = zeros(12,nDist);
    
    % Atomic position [Angstrom] of the magnetic atoms.
    atcoords       = obj.basisvector*atPos;
    basis_atcoords = obj.basisvector*mAtom.r;
    
    % Inside the original unit cell magn. atoms couple to atoms outside unit cell.
    index=1;
    for ii=1:nMagAtom
        for jj=1:nHalfCube
            sortM(12,index)  = norm(basis_atcoords(:,ii)-atcoords(:,jj));
            sortM(1:3,index) = atTr(:,jj);
            sortM(4,index)   = ii;
            sortM(5,index)   = atIndex(jj);
            index = index+1;
        end
    end
    
    % Couplings inside origin unit cell.
    for ii=1:(nMagAtom-1)
        for jj=(ii+1):nMagAtom
            sortM(12,index)   = norm(basis_atcoords(:,ii)-basis_atcoords(:,jj));
            sortM(1:3,index) = [0 0 0];
            sortM(4,index)   = ii;
            sortM(5,index)   = jj;
            index = index+1;
        end
    end
    % sort according to distance and cut away large distances
    sortM = sortrows(sortM(:,sortM(12,:)<param.maxDistance)',12)';
    
    % Finds the equivalent distances and index them in coupling.idx
    sortM(6,:) = cumsum([1 (sortM(12,2:end)-sortM(12,1:(end-1))) > tol]);
    
    aniso = int32(zeros(1,nMagAtom));
    % symmetry equivalent couplings
    if param.sym
        % get the symmetry operators
        [symOp, ~] = sw_gencoord(obj.lattice.sym);
        % store the final sorted couoplings in newM
        newM = zeros(6,0);
        ii  = 1;
        idx = 1;
        while (ii < max(sortM(6,:)))
            % select columns from sorM with a certain idx value
            sortMs = sortM(:,sortM(6,:) == ii);
            while (size(sortMs,2)>0)
                % atom positions and translation vector of the first atom
                r1 = mAtom.r(:,sortMs(4,1));
                r2 = mAtom.r(:,sortMs(5,1));
                dl = sortMs(1:3,1);
                % generate new atomic positions and translation vectors
                r1new = permute(mmat(symOp,r1),[1 3 2]);
                r2new = permute(mmat(symOp,r2),[1 3 2]);
                dlnew = permute(mmat(symOp,dl),[1 3 2]);
                % determine the new indices in mAtom
                [~, atom1] = isnew(mAtom.r,r1new,param.tol);
                [~, atom2] = isnew(mAtom.r,r2new,param.tol);
                % remove identical couplings
                
                
                % generate coupling columns
                genC = [[dlnew; atom1;atom2] [-dlnew;atom2;atom1]];
                % detemine the independent couplings
                symIdx = isnew2(genC,sortMs(1:5,:),param.tol);
                % move the non-unique (not new) couplings (symmetry equivalent ones)
                newM = [newM [sortMs(1:5,symIdx);ones(1,numel(symIdx))*idx]]; %#ok<AGROW>
                idx  = idx + 1;
                % remove the symmetry equivalent couplings
                sortMs(:,symIdx) = [];
            end
            ii = ii + 1;
        end
    else
        newM = sortM;
    end
else
    % If there is no magnetic atom the coupling and anisotropy are empty.
    newM  = zeros(12,0);
    aniso = int32(zeros(1,0));
    
end

coupling.dl      = int32(newM(1:3,:));
coupling.atom1   = int32(newM(4,:));
coupling.atom2   = int32(newM(5,:));
coupling.idx     = int32(newM(6,:));
coupling.mat_idx = int32(zeros(3,size(coupling.idx,2)));

obj.coupling         = coupling;
obj.single_ion.aniso = aniso;

validate(obj);

end

function [isnew, symIdx] = isnew(A,B, tol)
% selects the new vectors from B. Dimensions of A and B have to be [3 nA]
% and [3 nB] respectively.
% A vector in B is considered new, if d(mod(vA-vB,1))<tol.
%
% symIdx    Indices to re

nA = size(A,2);
nB = size(B,2);

notequal = sum(mod(abs(repmat(permute(A,[2 3 1]),[1 nB 1]) - repmat(permute(B,[3 2 1]),[nA 1 1])),1).^2,3)>tol.^2;

isnew = all(notequal,1);

idx = 1:nB;

symIdx = arrayfun(@(idx)find(~notequal(:,idx),1,'first'),idx(~isnew));

end

function symIdx = isnew2(A,B, tol)
% selects the new vectors from B. Dimensions of A and B have to be [3 nA]
% and [3 nB] respectively.
% A vector in B is considered new, if d(mod(vA-vB,1))<tol.
%
% symIdx    Indices to re

nA = size(A,2);
nB = size(B,2);

notequal = sum(abs(repmat(permute(A,[2 3 1]),[1 nB 1]) - repmat(permute(B,[3 2 1]),[nA 1 1])).^2,3)>tol.^2;

isnew = all(notequal,2);

idx = 1:nA;

symIdx = arrayfun(@(idx)find(~notequal(idx,:),1,'first'),idx(~isnew));

end

