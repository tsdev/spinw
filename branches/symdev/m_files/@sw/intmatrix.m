function [SS, SI, RR] = intmatrix(obj, varargin)
% creates the interactions matrices (connectors and values)
%
% [SS, SI, RR] = INTMATRIX(obj)
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
% For speedup (if atomic positions does not change) call 
% obj.intmatrix(true).
%

% Create parameters of magnetic atoms in the unit cell.
mAtom     = obj.matom(varargin{:});
nMagAtom  = size(mAtom.r,2);

coupling  = obj.coupling;

SS.all = double([coupling.dl; coupling.atom1; coupling.atom2]);

% Remove elements of coupling where mat_idx == 0.
Col1 = find(coupling.mat_idx(1,:) ~= 0);
Col2 = find(coupling.mat_idx(2,:) ~= 0);
Col3 = find(coupling.mat_idx(3,:) ~= 0);

JJ.idx = [coupling.mat_idx(1,Col1) coupling.mat_idx(2,Col2) coupling.mat_idx(3,Col3)];

SS.all = SS.all(:,[Col1 Col2 Col3]);

JJ.mat  = obj.matrix.mat(:,:,JJ.idx);

% For non P1 symmetry, calculate the interaction matrices
% TODO





% don't calculate these for speedup in case of fitting
if (nargin==1) || (nargin>1 && ~varargin{1})
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
SS.all = SS.all(:,sum(SS.all(6:end,:).^2,1)~=0);

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