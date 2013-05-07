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
%               the couplings for plotting.
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
%               spin(:)*A(:,:,i)*spin(:)'.
% SI.field      External magnetic field [Bx By Bz].
%
% RR            Positions of the atoms in lattice units, dimensions are
%               [3 nMAgExt].
%

inpForm.fname  = {'fitmode' 'plotmode'};
inpForm.defval = {0          false    };
inpForm.size   = {[1 1]      [1 1]    };

param = sw_readparam(inpForm, varargin{:});

% Create parameters of magnetic atoms in the unit cell.
mAtom    = obj.matom(param.fitmode);
nMagAtom = size(mAtom.r,2);
mat      = obj.matrix.mat;
nMat     = size(mat,3);

coupling = obj.coupling;

SS.all = double([coupling.dl; coupling.atom1; coupling.atom2; coupling.idx]);

if param.fitmode > 0
    % sum up couplings on the same bond
    
    mat = cat(3,mat,zeros(3));
    % Remove couplings where all mat_idx == 0.
    colSel = any(coupling.mat_idx ~= 0,1);
    
    JJ.idx = coupling.mat_idx(:,colSel);
    JJ.idx(JJ.idx == 0) = nMat + 1;
    SS.all = SS.all(:,colSel);
    SS.all(end+1,:) = 0;
    
    % sum the interactions on the same coupling
    JJ.mat  = mat(:,:,JJ.idx(1,:)) + mat(:,:,JJ.idx(2,:)) + mat(:,:,JJ.idx(3,:));
    
else
    % just keep all non-zero coupling
    mat_idx = coupling.mat_idx';
    JJ.idx  = mat_idx(mat_idx(:) ~= 0);
    
    colSel  = [find(coupling.mat_idx(1,:)~=0)' find(coupling.mat_idx(2,:)~=0)' find(coupling.mat_idx(3,:)~=0)'];
    SS.all = SS.all(:,colSel(:));
    
    SS.all = [SS.all; double(JJ.idx')];
    
    % sum the interactions on the same coupling
    JJ.mat  = mat(:,:,JJ.idx);
    
end

% For non P1 symmetry, calculate the interaction matrices
if obj.lattice.sym > 1
    % first positions of the couplings with identical idx values used to
    % generate the coupling matrices for the rest
    firstC = SS.all(1:5,[true logical(diff(SS.all(6,:)+100*SS.all(7,:)))]);
    
    % center positions of the first element of couplings with identical idx
    % values
    cPos = (mAtom.r(:,firstC(4,:)) + mAtom.r(:,firstC(5,:)) + firstC(1:3,:))/2;
    % operators that rotate the coupling matrices
    [~,~,~,rotOp] = sw_genatpos(obj.lattice.sym,cPos);
    
    % convert rotation operators to xyz Cartesian coordinate system
    A = obj.basisvector;
    rotOp = mmat(A,mmat(rotOp,inv(A)));
    
    % rotate the matrices: R*M*R'
    JJ.mat = mmat(rotOp,mmat(JJ.mat,permute(rotOp,[2 1 3])));
    
end

idxTemp = SS.all(6,:);
SS.all  = SS.all(1:5,:);

% don't calculate these for speedup in case of fitting
if param.fitmode > 1
    JJ.type = sw_mattype(JJ.mat);
    
    % Select only the isotropic exchange. Remove zero value elements.
    JJ.iso = squeeze(JJ.mat(1,1,:))';
    SS.iso = [SS.all(:,JJ.type == 1); JJ.iso(1,JJ.type == 1)];
    SS.iso = SS.iso(:,SS.iso(6,:)~=0);
    
    % Select only the anisotropic exchange. Remove zero value elements.
    JJ.ani = [squeeze(JJ.mat(1,1,:))'; squeeze(JJ.mat(2,2,:))'; squeeze(JJ.mat(3,3,:))'];
    SS.ani = [SS.all(:,JJ.type == 2); JJ.ani(:,JJ.type == 2)];
    SS.ani = SS.ani(:,sum(SS.ani(6:end,:).^2,1)~=0);
    
    % Select DM interactions. Remove zero value elements.
    JJ.dm = [squeeze(JJ.mat(2,3,:))'; squeeze(JJ.mat(3,1,:))'; squeeze(JJ.mat(1,2,:))'];
    SS.dm = [SS.all(:,JJ.type == 3); JJ.dm(:,JJ.type == 3)];
    SS.dm = SS.dm(:,sum(SS.dm(6:end,:).^2,1)~=0);
    
    % Select general interactions. Remove zero value elements.
    JJ.gen = reshape(JJ.mat,1,9,size(JJ.mat,3));
    SS.gen = [SS.all(:,JJ.type == 4); shiftdim(JJ.gen(1,:,JJ.type == 4),1)];
    SS.gen = SS.gen(:,sum(SS.gen(6:end,:).^2,1)~=0);
end

% Saves the whole interaction matrices into SS.all
SS.all = [SS.all; shiftdim(reshape(JJ.mat,1,9,size(JJ.mat,3)),1)];
SS.all = SS.all(:,sum(SS.all(6:end,:).^2,1) > 1e-10);

if param.plotmode
    % Saves all coupling matrix indices in SS.all in case of non-fitting mode
    % in the bottom row
    SS.all = [SS.all; double(JJ.idx'); idxTemp];
end

% Anisotropy matrix.
SI.aniso = zeros(3,3,nMagAtom);

if size(obj.single_ion.aniso,2) == nMagAtom
    for ii = 1:nMagAtom
        idx = obj.single_ion.aniso(ii);
        if any(idx)
            SI.aniso(:,:,ii) = obj.matrix.mat(:,:,idx);
            SI.aniso(:,:,ii) = (SI.aniso(:,:,ii)+SI.aniso(:,:,ii)')/2;
        end
    end
end

% Extend the lattice.
[mAtom, SI.aniso, SS] = sw_extendlattice(double(obj.mag_str.N_ext), mAtom, SI.aniso, SS);

% Save external field.
SI.field = obj.single_ion.field;

% Save the position of all atoms
RR = mAtom.RRext;

end