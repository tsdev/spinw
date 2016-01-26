function [SS, SI, RR] = intmatrix(obj, varargin)
% creates the interactions matrices (connectors and values)
%
% [SS, SI, RR] = INTMATRIX(obj, 'Option1', Value1, ...)
%
% Input:
%
% obj           Input sw class object.
%
% Options:
%
% fitmode       Can be used to speed up calculation, modes:
%               0   No speedup, default.
%               1   Only atomic positions are precalculated and equivalent
%                   coupling matrices are summed up.
%               2   Same as mode == 1, moreover only SS.all is calculated.
%               3   Same as mode == 2, but save symmetry information to sw
%                   object to reuse on next call with fitmode 3 (if other
%                   fitmode, or no fitmode used, recalculate this info).
% plotmode      If true, additional rows are added to SS.all, to identify
%               the couplings for plotting. Default is false.
% sortDM        If true each coupling is sorted for consistent plotting of
%               the DM interaction. Sorting is based on the dR distance
%               vector, pointing from atom1 to atom2. Its components should
%               fulfill the following rules in hierarchical order:
%                   1. dR(x) > 0
%                   2. dR(y) > 0
%                   3. dR(z) > 0.
%               Default is false.
% zeroC         Whether to output bonds with assigned matrices that are
%               zero. Default is false.
% extend        If true, all bonds in the magnetic supercell will be
%               generated, if false, only the bonds in the crystallographic
%               unit cell is calculated. Default is true.
% conjugate     Introduce the conjugate of the couplings (atom1 and atom2
%               exchanged). Default is false.
% rotMat        Rotate the J and A matrices according to the point group
%               operations between symmetry equivalent sites. Default is
%               true.
%
% Output:
%
% SS            Structure with  fields {iso,aniso,dm,gen,bq,dip}. It
%               describes the interactions between spins. Every field is a
%               matrix, where every column is a coupling between two spins.
%               The first 3 rows contain the unit cell translation vector
%               between the interacting spins, the 4th and 5th row contains
%               the indices of the two interacting spins in the 'spin'
%               variable. The following rows contains the strength of the
%               interaction. For isotropic exchange it is a single number,
%               for DM interaction [DMx; DMy; DMz], for anisotropic
%               interaction [Jxx; Jyy; Jzz] and for general interaction
%               [Jxx; Jxy; Jxz; Jyx; Jyy; Jyz; Jzx; Jzy; Jzz] and for
%               biquadratic exchange it is also a single number.
%               For example:
%                SS.iso = [dLatX; dLatY; dLatZ; spinIdx1; spinIdx2; Jval].
%               For plotmode true, two additional rows are added to SS.all,
%               that contains the idx indices of the obj.matrix(:,:,idx)
%               corresponding matrix for each coupling and the .idx values
%               of the couplings. The dip field contains the dipolar
%               interactions that are not added to the SS.all field.
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

%if obj.symbolic && obj.symmetry
%    if any(sw_mattype(obj.matrix.mat)~=1)
%        warning('sw:intmatrix:symmetry',['The non-isotropic symbolic matrices '...
%            'will be rotated unsing the point group operators, the result can be ugly!']);
%    end
%end

inpForm.fname  = {'fitmode' 'plotmode' 'zeroC' 'extend' 'conjugate' 'rotMat' 'sortDM'};
inpForm.defval = {0          false     false   true     false       true     false   };
inpForm.size   = {[1 1]      [1 1]     [1 1]   [1 1]    [1 1]       [1 1]    [1 1]   };

param = sw_readparam(inpForm, varargin{:});

% Create parameters of magnetic atoms in the unit cell.
mAtom    = obj.matom;
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
SS.all   = double([coupling.dl; coupling.atom1; coupling.atom2; coupling.idx]);

if param.fitmode > 0
    % sum up couplings on the same bond
    
    mat_type = double(coupling.type);
    
    % loop over the different types of couplings
    typeList = unique(mat_type(:));
    
    JJ.mat = zeros(3,3,0);
    JJ.idx = zeros(1,0);
    all_temp = zeros(8,0);
    
    for ii = 1:numel(typeList)
        mat_idx  = coupling.mat_idx;
        
        if numel(typeList)>1
            mat_idx(mat_type~=typeList(ii)) = 0;
        end
        % Remove couplings where all mat_idx == 0.
        colSel = any(mat_idx ~= 0,1);
        
        mat_temp = mat_idx(:,colSel);
        mat_temp(mat_temp == 0) = nMat + 1;
        % add 0 as last row instead of JJ.idx, since we sum up the matrices
        all_temp = [all_temp [SS.all(:,colSel);repmat([0; typeList(ii)],1,sum(colSel))]]; %#ok<AGROW>
        
        % sum the interactions on the same coupling
        JJ.mat = cat(3,JJ.mat,mat(:,:,mat_temp(1,:)) + mat(:,:,mat_temp(2,:)) + mat(:,:,mat_temp(3,:)));
        JJ.idx = [JJ.idx; mat_idx(mat_idx(:) ~= 0)];
    end
    SS.all = all_temp;
else
    % just keep all coupling that has an assigned matrix
    % TODO
    mat_idx  = coupling.mat_idx';
    mat_type = double(coupling.type)';
    
    JJ.idx  = mat_idx(mat_idx(:) ~= 0);
    
    colSel  = [find(coupling.mat_idx(1,:)~=0) find(coupling.mat_idx(2,:)~=0) find(coupling.mat_idx(3,:)~=0)];
    SS.all = SS.all(:,colSel(:));
    
    SS.all = [SS.all; double(JJ.idx')];
    
    % select the non-zero interactions into JJ.mat
    JJ.mat  = mat(:,:,JJ.idx);
    
    % add an extra type row to SS.all
    SS.all = [SS.all; mat_type(mat_idx(:) ~= 0)'];
end

% For non P1 symmetry and when the couplings are generated using symmetry,
% generate the Hamiltonian using the symmetry operators
if obj.sym
    % Generate the transformation matrix between lattice units and xyz
    % coordinate system
    A = obj.basisvector(false,obj.symbolic);
    
    % Generate anisotropy matrice using the space group symmetry
    if obj.symb
        nzeroA = ~sw_always(sum(SI.aniso(:).^2,1)==0);
    else
        nzeroA = any(SI.aniso(:));
    end
    
    % generate the rotation matrices
    if param.fitmode == 3
        if ~isfield(obj.lattice,'rotOp')
            [~, ~, ~, rotOp] = sw_genatpos(obj.lattice.sym,obj.unit_cell.r(:,~sw_always(obj.unit_cell.S==0)));
            obj.lattice.rotOp = rotOp;
        end
        rotOp = obj.lattice.rotOp;
    else
        if obj.symbolic
            [~, ~, ~, rotOp] = sw_genatpos(obj.lattice.sym,obj.unit_cell.r(:,~sw_always(obj.unit_cell.S==0)));
            % try to convert the rotation matrices into symbolic (only works
            % for matrices with 1s and 0s)
            %rotOp = round(rotOp);
            obj.lattice.rotOp = rotOp;
        else
            [~, ~, ~, rotOp] = sw_genatpos(obj.lattice.sym,obj.unit_cell.r(:,obj.unit_cell.S>0));
        end
        if isfield(obj.lattice,'rotOp')
            rmfield(obj.lattice,'rotOp');
        end
    end
    % convert rotation operators to xyz Cartesian coordinate system
    %if sum(abs((A-diag(diag(A)))/max(A))) > 1e-15
    rotOp = mmat(A,mmat(rotOp,inv(A)));
    %end
    
    if param.rotMat
        if nzeroA
            % rotate the matrices: R*M*R'
            SI.aniso = mmat(rotOp,mmat(SI.aniso,permute(rotOp,[2 1 3])));
        end
        % Generate g-tensor using the space group symmetry
        % rotate the matrices: R*M*R'
        SI.g = mmat(rotOp,mmat(SI.g,permute(rotOp,[2 1 3])));
    end
    
    % Generate interaction matrices using the space group symmetry
    if ~isempty(SS.all)
        % first positions of the couplings with identical idx values used to
        % generate the coupling matrices for the rest
        firstC = SS.all(1:5,[true logical(diff(SS.all(6,:)+100*SS.all(7,:)))]);
        % produce the space group symmetry operators
        if param.fitmode == 3
            if isfield(obj.lattice,'symOp')
                rotOp = obj.lattice.symOp;
            else
                [symOp, symTr] = sw_gencoord(obj.lattice.sym);
                % select rotation matrices for each generated coupling
                rotOp = zeros(3,3,0);
                for ii = 1:size(firstC,2)
                    [~, rotIdx] = sw_gensymcoupling(obj, firstC(:,ii), {symOp, symTr}, 1e-5);
                    rotOp = cat(3,rotOp,symOp(:,:,rotIdx));
                end
        
                %if obj.symbolic
                %    rotOp = round(rotOp);
                %end
        
                % convert rotation operators to xyz Cartesian coordinate system
                %if sum(abs((A-diag(diag(A)))/max(A))) > 1e-15
                rotOp = mmat(A,mmat(rotOp,inv(A)));
                %end
                obj.lattice.symOp = rotOp;
            end
        else
            [symOp, symTr] = sw_gencoord(obj.lattice.sym);
            % select rotation matrices for each generated coupling
            rotOp = zeros(3,3,0);
            for ii = 1:size(firstC,2)
                [~, rotIdx] = sw_gensymcoupling(obj, firstC(:,ii), {symOp, symTr}, 1e-5, true);
                rotOp = cat(3,rotOp,symOp(:,:,rotIdx));
            end
            rotOp = mmat(A,mmat(rotOp,inv(A)));
            if isfield(obj.lattice,'symOp') 
                rmfield(obj.lattice,'symOp');
            end
        end
        % rotate the matrices: R*M*R'
        if param.rotMat
            JJ.mat = mmat(rotOp,mmat(JJ.mat,permute(rotOp,[2 1 3])));
        end
    end
end

mat_type = SS.all(8,:);
idxTemp  = SS.all(6,:);
SS.all   = SS.all(1:5,:);

% don't calculate these for speedup in case of fitting
if param.fitmode < 2
    JJ.type = sw_mattype(JJ.mat);
    
    % new type for biquadratic exchange
    if any(JJ.type(mat_type==1)~=1)
        error('sw:intmatrix:DataError','Biquadratic exchange matrix has to be isotropic!')
    end
    % for biquadratic exchange type = 5
    JJ.type(mat_type==1) = 5;
    
    % Select only the isotropic exchange. Remove zero value elements.
    JJ.iso = squeeze(JJ.mat(1,1,:))';
    SS.iso = [SS.all(:,JJ.type == 1); JJ.iso(1,JJ.type == 1)];
    
    if ~isempty(SS.iso)
        SS.iso = SS.iso(:,SS.iso(6,:)~=0);
    end
    
    % Select only the isotropic exchange. Remove zero value elements.
    JJ.bq = squeeze(JJ.mat(1,1,:))';
    SS.bq = [SS.all(:,JJ.type == 5); JJ.iso(1,JJ.type == 5)];
    
    if ~isempty(SS.bq)
        SS.bq = SS.bq(:,SS.bq(6,:)~=0);
    end
    
    % Select only the anisotropic exchange. Remove zero value elements.
    JJ.ani = [squeeze(JJ.mat(1,1,:))'; squeeze(JJ.mat(2,2,:))'; squeeze(JJ.mat(3,3,:))'];
    SS.ani = [SS.all(:,JJ.type == 2); JJ.ani(:,JJ.type == 2)];
    
    if obj.symbolic
        idx = any(~sw_always(SS.ani(6:8,:)==0));
    else
        idx = any(SS.ani(6:8,:));
    end
    SS.ani = SS.ani(:,idx);
    
    % Select DM interactions. Remove zero value elements.
    JJ.dm = [squeeze(JJ.mat(2,3,:))'; squeeze(JJ.mat(3,1,:))'; squeeze(JJ.mat(1,2,:))'];
    SS.dm = [SS.all(:,JJ.type == 3); JJ.dm(:,JJ.type == 3)];
    
    if obj.symbolic
        idx = any(~sw_always(SS.dm(6:8,:)==0));
    else
        idx = any(SS.dm(6:8,:));
    end
    
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
SS.all = [SS.all; shiftdim(reshape(JJ.mat,1,9,size(JJ.mat,3)),1);mat_type];

if obj.coupling.rdip > 0
    % add dipolar interactions
    
    % calculate bond vectors in xyz coordinate system
    dr   = double(coupling.dl)+mAtom.r(:,coupling.atom2)-mAtom.r(:,coupling.atom1);
    drA  = obj.basisvector*dr;
    lA   = sqrt(sum(drA.^2,1));
    % calculate the matrices only for bond vector shorter than rdip
    rSel = lA<obj.coupling.rdip;
    nR   = sum(rSel);
    % remove longer bonds
    rAN  = bsxfun(@rdivide,drA(:,rSel),lA(1,rSel));
    
    rrmat = bsxfun(@times,permute(rAN,[1 3 2]),permute(rAN,[3 1 2])) - repmat(eye(3),[1 1 nR]);
    
    Edip = obj.unit.mu0*obj.unit.muB^2/4/pi;
    % dipole-dipole interaction matrices
    rrmat = bsxfun(@times,permute(Edip./lA(1,rSel).^3,[1 3 2]),rrmat);
    % multiply with the g-tensors on the left and right sides
    rrmat = mmat(permute(SI.g(:,:,coupling.atom1(rSel)),[2 1 3]),rrmat);
    Jdip  = mmat(rrmat,SI.g(:,:,coupling.atom2(rSel)));
    
    % create the matrix for dipolar interaction list
    SS.dip = double([coupling.dl(:,rSel); coupling.atom1(1,rSel); ...
        coupling.atom2(1,rSel)]);
    SS.dip = [SS.dip; reshape(Jdip,9,[]); zeros(1,nR)];
else
    SS.dip = zeros(15,0);
end

% Cut out zero valued interactions
if ~param.zeroC
    if obj.symb
        nzeroJ = ~sw_always(sum(SS.all(6:14,:).^2,1)==0);
    else
        nzeroJ  = sum(SS.all(6:14,:).^2,1) > 1e-10;
    end
    SS.all  = SS.all(:,nzeroJ);
    JJ.idx  = JJ.idx(nzeroJ);
    JJ.mat  = JJ.mat(:,:,nzeroJ);
    idxTemp = idxTemp(nzeroJ);
end

if param.plotmode
    % Saves all coupling matrix indices in SS.all in case of non-fitting mode
    % in the bottom row
    SS.all   = [SS.all(1:14,:); double(JJ.idx'); idxTemp;SS.all(15,:)];
end

if param.sortDM && (~isempty(SS.all))
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

if param.extend
    % Extend the lattice for magnetic interactions
    nExt = double(obj.mag_str.N_ext);
    [mAtom, SS] = sw_extendlattice(nExt, mAtom, SS);
    SI.aniso = repmat(SI.aniso, [1 1 prod(nExt)]);
    SI.g     = repmat(SI.g, [1 1 prod(nExt)]);
    
    % Save the position of all atoms
    RR = mAtom.RRext;
else
    RR = mAtom.r;
end

if param.conjugate
    % Introduce the opposite couplings.
    % (i-->j) and (j-->i)
    % transpose the JJ matrix as well [1 2 3 4 5 6 7 8 9] --> [6 9 12 7 10 13 8 11 14]
    % this step is not necessary for diagonal exchange matrices and
    % biquadratic exchange
    if numel(SS.all) > 0
        new         = [SS.all(1:3,:)   -SS.all(1:3,:)  ];
        new(4:5,:)  = [SS.all([4 5],:)  SS.all([5 4],:)];
        new(6:14,:) = [SS.all(6:14,:)   SS.all([6 9 12 7 10 13 8 11 14],:)]/2;
        new(15,:)   = [SS.all(end,:)    SS.all(end,:)];
        SS.all      = new;
    end
    
    if numel(SS.dip) > 0
        new         = [SS.dip(1:3,:)   -SS.dip(1:3,:)  ];
        new(4:5,:)  = [SS.dip([4 5],:)  SS.dip([5 4],:)];
        new(6:14,:) = [SS.dip(6:14,:)   SS.dip([6 9 12 7 10 13 8 11 14],:)]/2;
        new(15,:)   = [SS.dip(end,:)    SS.dip(end,:)];
        SS.dip      = new;
    end
end

% Save external field.
SI.field = obj.single_ion.field;

end
