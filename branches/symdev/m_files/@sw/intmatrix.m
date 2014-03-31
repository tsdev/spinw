function [SS, SI, RR] = intmatrix(obj, varargin)
% creates the interactions matrices (connectors and values)
%
% [SS, SI, RR] = INTMATRIX(obj, 'Option1', Value1, ...)
%
% Options:
%
% fitmode       To speed up calculation, modes:
%               1   only atomic positions are precalculated and equivalent
%                   coupling matrices are summed up
%               2   as mode == 1, moreover only SS.all is calculated.
% plotmode      If true, additional rows are added to SS.all, to identify
%               the couplings for plotting and each coupling is sorted for
%               consistent plotting of the DM interaction. Sorting is based
%               on the dr distance vector, pointing from atom1 to atom2.
%               Its components should fulfill the following rules in
%               hierarchical order:
%                   1. dr(x) > 0
%                   2. dr(y) > 0
%                   3. dr(z) > 0.
% zeroC         Whether to give couplings with assigned matrices that are
%               zero. Default is false.
% extend        If true, all bonds in the magnetic supercell will be
%               generated, if false, only the bonds in the crystallographic
%               unit cell is calculated. Default is true.
%
% Output:
%
% obj           Input onject contains structural data, sw type.
% SS            Structure with  fields {iso,aniso,dm,gen}. It describes
%               the interactions between spins. Every field is a matrix,
%               where every column is a coupling between two spins. The
%               first 3 rows contain the unit cell translation vector
%               between the interacting spins, the 4th and 5th row contains
%               the indices of the two interacting spins in the 'spin'
%               variable. The following rows contains the strength of the
%               interaction. For isotropic exchange it is a single number,
%               for DM interaction [DMx; DMy; DMz], for anisotropic
%               interaction [Jxx; Jyy; Jzz] and for general interaction
%               [Jxx; Jxy; Jxz; Jyx; Jyy; Jyz; Jzx; Jzy; Jzz]
%               For example:
%                SS.iso = [dLatX; dLatY; dLatZ; spinIdx1; spinIdx2; Jval].
%               For plotmode true, two additional rows are added to SS.all,
%               that contains the idx indices of the obj.matrix(:,:,idx)
%               corresponding matrix for each coupling and the .idx values
%               of the couplings.
%
% SI            Single ion energy, due to anisotropy and magnetic field.
% SI.aniso      Matrix with dimensions of [3 3 nMagAtom] sized matrix,
%               where the energy of the i-th spin is
%               E_aniso = spin(:)*A(:,:,i)*spin(:)'.
% SI.g          g-tensor, with dimensions of [3 3 nMagAtom]. It determines
%               the energy of the magnetic moment in external field:
%               E_field = B(:)*g(:,:,i)*spin(:)'.
% SI.field      External magnetic field [Bx By Bz].
%
% RR            Positions of the atoms in lattice units, dimensions are
%               [3 nMAgExt].
%
% See also SW.COUPLINGTABLE.
%

inpForm.fname  = {'fitmode' 'plotmode' 'zeroC' 'extend'};
inpForm.defval = {0          false     false   true    };
inpForm.size   = {[1 1]      [1 1]     [1 1]   [1 1]   };

param = sw_readparam(inpForm, varargin{:});

% Create parameters of magnetic atoms in the unit cell.
mAtom    = obj.matom(param.fitmode);
nMagAtom = size(mAtom.r,2);
mat      = obj.matrix.mat;
nMat     = size(mat,3);

% Add extra zero matrix to the end of the matrix list
mat = cat(3,mat,zeros(3));
% Add another extra matrix for g=2 default tensor
mat = cat(3,mat,2*eye(3));

% Anisotropy matrix.
if size(obj.single_ion.aniso,2) == nMagAtom
    idx = obj.single_ion.aniso;
    % the non-assigned atoms get a pointer to the last zero energy matrix
    idx(idx == 0) = nMat + 1;
    SI.aniso = mat(:,:,idx);
    % keeping only the symmetric part of the anisotropy
    SI.aniso = (SI.aniso + permute(SI.aniso,[2 1 3]))/2;
else
    SI.aniso = zeros(3,3,nMagAtom);
end

% g-tensor.
if size(obj.single_ion.g,2) == nMagAtom
    idx = obj.single_ion.g;
    % the non-assigned atoms get a pointer to the last zero energy matrix
    idx(idx == 0) = nMat + 2;
    SI.g = mat(:,:,idx);
    % keeping only the symmetric part of the g-tensor
    SI.g = (SI.g + permute(SI.g,[2 1 3]))/2;
else
    % default g-tensor value of 2
    SI.g = repmat(2*eye(3),[1 1 nMagAtom]);
end

% Couplings
coupling = obj.coupling;
SS.all = double([coupling.dl; coupling.atom1; coupling.atom2; coupling.idx]);

% just keep all non-zero coupling
mat_idx = coupling.mat_idx;
    
if param.fitmode > 0
    % sum up couplings on the same bond
    
    % Remove couplings where all mat_idx == 0.
    colSel = any(mat_idx ~= 0,1);
    
    JJ.idx = coupling.mat_idx(:,colSel);
    JJ.idx(JJ.idx == 0) = nMat + 1;
    SS.all = SS.all(:,colSel);
    SS.all(end+1,:) = 0;
    
    % sum the interactions on the same coupling
    JJ.mat  = mat(:,:,JJ.idx(1,:)) + mat(:,:,JJ.idx(2,:)) + mat(:,:,JJ.idx(3,:));
    JJ.idx  = mat_idx(mat_idx(:) ~= 0);
else
    % just keep all non-zero coupling
    %mat_idx = coupling.mat_idx';
    JJ.idx  = mat_idx(mat_idx(:) ~= 0);
    
    colSel  = [find(coupling.mat_idx(1,:)~=0) find(coupling.mat_idx(2,:)~=0) find(coupling.mat_idx(3,:)~=0)];
    SS.all = SS.all(:,colSel(:));
    
    SS.all = [SS.all; double(JJ.idx(:)')];
    
    % select the non-zero interactions into JJ.mat
    JJ.mat  = mat(:,:,JJ.idx);
    
end

% For non P1 symmetry and when the couplings are generated using symmetry,
% generate the Hamiltonian using the symmetry operators
if obj.sym
    % Generate the transformation matrix between lattice units and xyz
    % coordinate system
    A = obj.basisvector;
    
    % Generate anisotropy matrice using the space group symmetry
    if obj.symb
        nzeroA = sum(SI.aniso(:).^2,1)==0;
        if ~isa(nzeroA,'logical')
            nzeroA = isAlways(nzeroA);
        end
        nzeroA = ~nzeroA;
    else
        nzeroA  = any(SI.aniso(:));
    end
    
    if nzeroA
        [~, ~, ~, rotOp] = sw_genatpos(obj.lattice.sym,obj.unit_cell.r(:,obj.unit_cell.S>0));
        % convert rotation operators to xyz Cartesian coordinate system
        rotOp = mmat(A,mmat(rotOp,inv(A)));
        % rotate the matrices: R*M*R'
        SI.aniso = mmat(rotOp,mmat(SI.aniso,permute(rotOp,[2 1 3])));
    end
    
    % Generate g-tensor using the space group symmetry
    [~, ~, ~, rotOp] = sw_genatpos(obj.lattice.sym,obj.unit_cell.r(:,obj.unit_cell.S>0));
    % convert rotation operators to xyz Cartesian coordinate system
    rotOp = mmat(A,mmat(rotOp,inv(A)));
    % rotate the matrices: R*M*R'
    SI.g = mmat(rotOp,mmat(SI.g,permute(rotOp,[2 1 3])));
    
    
    % Generate interaction matrices using the space group symmetry
    if numel(SS.all) > 0
        % first positions of the couplings with identical idx values used to
        % generate the coupling matrices for the rest
        firstC = SS.all(1:5,[true logical(diff(SS.all(6,:)+100*SS.all(7,:)))]);
        % produce the space group symmetry operators
        [symOp, symTr] = sw_gencoord(obj.lattice.sym);
        rotOp = zeros(3,3,0);
        % select rotation matrices for each generated coupling
        for ii = 1:size(firstC,2)
            [~, rotIdx] = sw_gensymcoupling(obj, firstC(:,ii), {symOp, symTr}, 1e-5, true);
            rotOp = cat(3,rotOp,symOp(:,:,rotIdx));
        end
        % convert rotation operators to xyz Cartesian coordinate system
        rotOp = mmat(A,mmat(rotOp,inv(A)));
        
        % rotate the matrices: R*M*R'
        JJ.mat = mmat(rotOp,mmat(JJ.mat,permute(rotOp,[2 1 3])));
    end
end

idxTemp = SS.all(6,:);
SS.all  = SS.all(1:5,:);

% don't calculate these for speedup in case of fitting
if param.fitmode < 2
    JJ.type = sw_mattype(JJ.mat);
    
    % Select only the isotropic exchange. Remove zero value elements.
    JJ.iso = squeeze(JJ.mat(1,1,:))';
    SS.iso = [SS.all(:,JJ.type == 1); JJ.iso(1,JJ.type == 1)];

    SS.iso = SS.iso(:,SS.iso(6,:)~=0);
    
    % Select only the anisotropic exchange. Remove zero value elements.
    JJ.ani = [squeeze(JJ.mat(1,1,:))'; squeeze(JJ.mat(2,2,:))'; squeeze(JJ.mat(3,3,:))'];
    SS.ani = [SS.all(:,JJ.type == 2); JJ.ani(:,JJ.type == 2)];
    
    
    %idx = sum(SS.ani(6:end,:).^2,1)~=0;
    idx = any(SS.ani(6:end,:));
    SS.ani = SS.ani(:,idx);
    
    % Select DM interactions. Remove zero value elements.
    JJ.dm = [squeeze(JJ.mat(2,3,:))'; squeeze(JJ.mat(3,1,:))'; squeeze(JJ.mat(1,2,:))'];
    SS.dm = [SS.all(:,JJ.type == 3); JJ.dm(:,JJ.type == 3)];
    
    idx = any(SS.dm(6:end,:));
    %SS.dm = SS.dm(:,sum(SS.dm(6:end,:).^2,1)~=0);
    SS.dm = SS.dm(:,idx);
    
    % Select general interactions. Remove zero value elements.
    JJ.gen = reshape(JJ.mat,1,9,size(JJ.mat,3));
    SS.gen = [SS.all(:,JJ.type == 4); shiftdim(JJ.gen(1,:,JJ.type == 4),1)];
    
    idx = any(SS.gen(6:end,:));
    %SS.gen = SS.gen(:,sum(SS.gen(6:end,:).^2,1)~=0);
    SS.gen = SS.gen(:,idx);
    
end

% Saves the whole interaction matrices into SS.all
SS.all = [SS.all; shiftdim(reshape(JJ.mat,1,9,size(JJ.mat,3)),1)];

% Cut out zero valued interactions
if ~param.zeroC
    if obj.symb
        nzeroJ = sum(SS.all(6:end,:).^2,1)==0;
        if ~isa(nzeroJ,'logical')
            nzeroJ = isAlways(nzeroJ);
        end
        nzeroJ = ~nzeroJ;
    else
        nzeroJ  = sum(SS.all(6:end,:).^2,1) > 1e-10;
    end
    % for symbolic variables
    %nzeroJ  = isAlways(sum(SS.all(6:end,:).^2,1) > 1e-10);
    SS.all  = SS.all(:,nzeroJ);
    JJ.idx  = JJ.idx(nzeroJ);
    JJ.mat  = JJ.mat(:,:,nzeroJ);
    idxTemp = idxTemp(nzeroJ);
end

if param.plotmode
    % Saves all coupling matrix indices in SS.all in case of non-fitting mode
    % in the bottom row
    SS.all   = [SS.all; double(JJ.idx'); idxTemp];
    
    if ~isempty(SS.all)
        % sort properly the atom1-atom2 pairs for the DM interaction
        % based on the vector pointing from atom1 to atom2
        % rv = r_atom2 + dl - r_atom1
        rv = mAtom.r(:,SS.all(5,:)) + SS.all(1:3,:) - mAtom.r(:,SS.all(4,:));
        rmax = max(max(double(rv)));
        
        multL = ceil(fliplr(cumprod([1 [1 1]*(rmax+1)])));
        
        % find the couplings that have to be flipped
        flip = find(sum(bsxfunsym(@times,rv,multL'),1) < 0);
        % flip the selected couplings
        SS.all(1:3,flip)   = -SS.all(1:3,flip);
        SS.all([4 5],flip) =  SS.all([5 4],flip);
        % change the sign of the DM interaction (transpose J matrices)
        SS.all(6:14,flip)  = SS.all([1 4 7 2 5 8 3 6 9]+5,flip);
    end
    RR = mAtom.r;
end

if param.extend
    % Extend the lattice for magnetic interactions
    nExt = double(obj.mag_str.N_ext);
    [mAtom, SS] = sw_extendlattice(nExt, mAtom, SS);
    SI.aniso = repmat(SI.aniso, [1 1 prod(nExt)]);
    SI.g     = repmat(SI.g, [1 1 prod(nExt)]);
    
    % Save the position of all atoms
    RR = mAtom.RRext;
end

% Save external field.
SI.field = obj.single_ion.field;

end