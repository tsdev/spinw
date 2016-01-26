function gencoupling(obj, varargin)
% generates the COUPLING property of sw object
%
% GENCOUPLING(obj, 'option1', value, ...)
%
% The function calculates equivalent bonds between magnetic atoms. These
% are determined either based on crystal symmetry or bond length (with
% tolDist tolerance). If the space group index of 0 is defined
% (obj.lattice.sym=0), the equivalent bonds will be based on bond length.
% For space group index larger than 0, the symmetry equivalent bonds will
% be determined. This can ve overwritten by the forceNoSym parameter to
% consider bond length.
%
% IMPORTANT!
%   This function has to be used after the crystal structure is defined.
%   The SW.ADDCOUPLING, SW.COUPLINGTABLE functions will only work
%   afterwards.
%
% Input:
%
% obj           sw class object.
%
% Options:
%
% forceNoSym    If true, equivalent bonds are generated based on
%               bond length with .tolDist tolerance. If false symmetry
%               operators will be used if they are given
%               (obj.lattice.sym>0).
% nUnitCell     Edge length of the parallelepiped (same along a,b and c)
%               withing the algorithm searches for nearest neighbors in
%               lattice units. Default is 3.
% maxDistance   Maximum bond length that will be stored in the
%               obj.coupling property in units of Angstrom. Default is 8.0.
% maxBond       Maximum number of generated inequivalent bonds. If
%               zero, maxDistance is considered to limit the longest
%               bond. Default is zero.
% tolDist       Tolerance of distance, within two bonds are regarded
%               equivalent, default is 1e-3 Angstrom. Only used, when no
%               space group is defined.
% dMin          Minimum bond length, below which an error is triggered.
%               Default value is 0.5 Angstrom.
%
% Output:
%
% The obj.coupling field will be filled with values, depending on the
% crystal geometry.
%
% Example:
%
% cryst = sw;
% cryst.genlattice('lat_const',[3 3 5],'angled',[90 90 120])
% cryst.addatom('r',[0 0 0])
% cryst.gencoupling
% cryst.couplingtable(1:3)
%
% A triangular lattice is created in cryst and after using gencoupling()
% the couplingtable() function lists the 1st, 2nd and 3rd neighbor bonds.
%
% See also SW, SW.SYMMETRY, SW.NOSYM.
%

isSym = obj.lattice.sym > 0;

inpForm.fname  = {'forceNoSym' 'nUnitCell' 'maxDistance' 'tol' 'tolDist' 'dMin' 'maxBond'};
inpForm.defval = {false        3           8             1e-5  1e-3      0.5    0        };
inpForm.size   = {[1 1]        [1 1]       [1 1]         [1 1] [1 1]     [1 1]  [1 1]    };

param = sw_readparam(inpForm, varargin{:});
tol   = param.tol;
tolD  = param.tolDist;

% force no symmetry operator mode
if param.forceNoSym
    isSym = false;
end

% save the sym/nosym method into obj
obj.sym = isSym;

% Number of unit cells along any direction to find nearest neighbours.
nUnitCell = param.nUnitCell;

% Generate all atomic positions in the unit cell.
mAtom      = obj.matom;

% Number of magnetic atoms in the unit cell.
nMagAtom = size(mAtom.r,2);

if nMagAtom == 0
    error('sw:gencoupling:NoMagAtom','There is no magnetic atom in the unit cell with S>0!');
end


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
sortM = zeros(7,nDist);

% Atomic position [Angstrom] of the magnetic atoms.
atcoords       = obj.basisvector*atPos;
basis_atcoords = obj.basisvector*mAtom.r;

% Inside the original unit cell magn. atoms couple to atoms outside unit cell.
index=1;
for ii=1:nMagAtom
    for jj=1:nHalfCube
        sortM(7,index)  = norm(basis_atcoords(:,ii)-atcoords(:,jj));
        sortM(1:3,index) = atTr(:,jj);
        sortM(4,index)   = ii;
        sortM(5,index)   = atIndex(jj);
        index = index+1;
    end
end

% Couplings inside origin unit cell.
for ii=1:(nMagAtom-1)
    for jj=(ii+1):nMagAtom
        sortM(7,index)   = norm(basis_atcoords(:,ii)-basis_atcoords(:,jj));
        sortM(1:3,index) = [0 0 0];
        sortM(4,index)   = ii;
        sortM(5,index)   = jj;
        index = index+1;
    end
end
% sort according to distance and cut away large distances
sortM = sortrows(sortM(:,sortM(7,:)<=param.maxDistance)',7)';

% Finds the equivalent distances and index them in coupling.idx
sortM(6,:) = cumsum([1 (sortM(7,2:end)-sortM(7,1:(end-1))) > tolD]);

% check whether some atoms are too close
if sortM(7,1)<param.dMin
    error('sw:gencoupling:AtomPos',['Some atoms are too close (Dmin=' num2str(sortM(7,1)) '<' num2str(param.dMin) '), check your crystal structure!']);
end

% symmetry equivalent couplings
if isSym
    % get the symmetry operators
    [symOp, symTr] = sw_gencoord(obj.lattice.sym);
    % store the final sorted couoplings in newM
    newM = zeros(6,0);
    ii  = 1;
    idx = 1;
    while (ii <= max(sortM(6,:)))
        % select columns from sorM with a certain idx value
        sortMs = sortM(:,sortM(6,:) == ii);
        while (size(sortMs,2)>0)
            [genC, unC] = sw_gensymcoupling(obj, sortMs(:,1), {symOp, symTr}, tol);
            genCAll = [genC [-genC(1:3,:);genC([5 4],:)]];
            % remove from sortMs the identical couplings
            iNew = isnew(genCAll(1:5,:),sortMs(1:5,:),tol);
            sortMs(:,~iNew) = [];
            % remove identical couplings from the symmetry generated
            % list
            genC(:,~unC) = [];
            if sum(~iNew) ~= sum(unC)
                error('sw:gencoupling:SymProblem','Symmetry error! ii=%d, idx=%d. Try to change ''tol'' parameter.',ii,idx);
            end
            % move the non-unique (not new) couplings (symmetry equivalent ones)
            newM = [newM [genC;ones(1,size(genC,2))*idx]]; %#ok<AGROW>
            idx  = idx + 1;
        end
        ii = ii + 1;
    end
else
    newM = sortM;
end

% cut the number of bonds
if param.maxBond > 0
    newM = newM(:,newM(6,:)<=param.maxBond);
end

% default anisotropy and g-tensor values
aniso = int32(zeros(1,nMagAtom));
g     = int32(zeros(1,nMagAtom));

% save output structure
coupling.dl      = int32(newM(1:3,:));
coupling.atom1   = int32(newM(4,:));
coupling.atom2   = int32(newM(5,:));
coupling.idx     = int32(newM(6,:));
coupling.mat_idx = int32(zeros(3,size(coupling.idx,2)));
coupling.type    = int32(zeros(3,size(coupling.idx,2)));
coupling.rdip    = obj.coupling.rdip;

obj.coupling         = coupling;
obj.single_ion.aniso = aniso;
obj.single_ion.g     = g;

validate(obj);

end

function [isnew, symIdx] = isnew(A,B, tol)
% selects the new vectors from B. Dimensions of A and B have to be [3 nA]
% and [3 nB] respectively.
%
% A vector in B is considered new, if d(mod(vA-vB,1))<tol.
%

nA = size(A,2);
nB = size(B,2);

notequal = sum(abs(repmat(permute(A,[2 3 1]),[1 nB 1]) - repmat(permute(B,[3 2 1]),[nA 1 1])).^2,3) > tol;

isnew = all(notequal,1);

idx = 1:nB;

symIdx = arrayfun(@(idx)find(~notequal(:,idx),1,'first'),idx(~isnew));

end