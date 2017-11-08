function varargout = newcell(obj,varargin)
% transforms lattice
% 
% ### Syntax
% 
% `newcell(obj,Name,Value)`
%
% `T = newcell(obj,Name,Value)`
% 
% ### Description
% 
% `T = newcell(obj,Name,Value)` redefines the unit cell using new basis
% vectors. The input three basis vectors are in lattice units of the
% original cell and define a parallelepiped. The atoms from the original
% unit cell will fill the new unit cell and if the two cells are compatible
% the structure won't change. The magnetic structure, bonds and single ion
% property definitions will be erased. The new cell will have different
% reciprocal lattice, however the original reciprocal lattice units will be
% retained automatically. To use the new reciprocal lattice, set the
% `'keepq'` option to `false`. In the default case the [spinw.spinwave]
% function will calculate spin wave dispersion at reciprocal lattice points
% of the original lattice. The transformation between the two lattices is
% stored in `spinw.unit.qmat`.
% 
% ### Examples
% 
% In this example we generate the triangular lattice antiferromagnet and
% convert the hexagonal cell to orthorhombic. This doubles the number of
% magnetic atoms in the cell and changes the reciprocal lattice. However we
% set `'keepq'` parameter to `true` to able to index the reciprocal lattice
% of the orthorhombic cell with the reciprocal lattice of the original
% hexagonal cell. To show that the two models are equivalent, we calculate
% the spin wave spectrum on both model using the same rlu. On the
% orthorhombic cell, the $Q$ value will be converted automatically and the
% calculated spectrum will be the same for both cases.
%
% ```
% >>tri = sw_model('triAF',1)
% >>tri_orth = copy(tri)
% >>tri_orth.newcell('bvect',{[1 0 0] [1 2 0] [0 0 1]},'keepq',true)
% >>tri_orth.gencoupling
% >>tri_orth.addcoupling('bond',1,'mat','J_1')
% >>newk = ((tri_orth.unit.qmat)*tri.magstr.k')'
% >>tri_orth.genmagstr('mode','helical','k',newk,'S',[1 0 0]')
% >>plot(tri_orth)
% >>>swplot.zoom(1.5)
% >>snapnow
% >>>figure
% >>subplot(2,1,1)
% >>sw_plotspec(sw_egrid(tri.spinwave({[0 0 0] [1 1 0] 501})),'mode','color','dE',0.2)
% >>subplot(2,1,2)
% >>spec = tri_orth.spinwave({[0 0 0] [1 1 0] 501});
% >>sw_plotspec(sw_egrid(tri_orth.spinwave({[0 0 0] [1 1 0] 501})),'mode','color','dE',0.2)
% >>snapnow
% ```
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Name-Value Pair Arguments
% 
% `'bvect'`
% : Defines the new lattice vectors in the original lattice
%   coordinate system. Cell with the following elements
%   `{v1 v2 v3}` or a $[3\times 3]$ matrix with `v1`, `v2` and `v3` as column
%   vectors: `[v1 v2 v3]`. Default value is `eye(3)` for indentity
%   transformation.
% 
% `'bshift'`
% : Row vector that defines a shift of the position of the unit cell.
%   Default value is `[0 0 0]`.
% 
% `'keepq'`
% : If true, the reciprocal lattice units of the new model will be
%   the same as in the old model. This is achieved by storing the
%   transformation matrix between the new and the old coordinate system in
%   `spinw.unit.qmat` and applying it every time a reciprocal space
%   definition is invoked, such as in [spinw.spinwave]. Default value is
%   `false`.
% 
% ### Output Arguments
% 
% `T`
% : Transformation matrix that converts $Q$ points (in reciprocal
%       lattice units) from the old reciprocal lattice to the new
%       reciprocal lattice as follows:
%   ```
%   Qrlu_new = T * Qrlu_old
%   ```
%   where the $Q$ vectors are row vectors with 3 elements.
% 
% ### See Also
% 
% [spinw.genlattice] \| [spinw.gencoupling] \| [spinw.nosym]
%
% *[rlu]: reciprocal lattice unit
%

if nargin <= 1
    help spinw.newcell
    return
end

inpForm.fname  = {'bvect' 'bshift' 'keepq'};
inpForm.defval = {eye(3)  [0 0 0]  false  };
inpForm.size   = {[-1 3]  [1 3]    [1 1]  };

param = sw_readparam(inpForm, varargin{:});

if ~iscell(param.bvect) && size(param.bvect,1)~=3
    error('spinw:newcell:WrongInput','Input has to be 1x3 cell or 3x3 matrix!');
end

% shift
bshift = param.bshift(:);

% here 3 coordinate systems are used:
% - xyz real space, Cartesian
% - original unit cell a,b and c vectors
% - new unit cell a',b' and c' vectors


% transformation matrix from the new lattice units into the original
% lattice
if iscell(param.bvect)
    Tn_o = [param.bvect{1}(:) param.bvect{2}(:) param.bvect{3}(:)];
else
    Tn_o = param.bvect;
end

% transformation from the original lattice into xyz real space
% xyz = To_xyz * v_orig
To_xyz = obj.basisvector;

% the new basis vectors
basisvector2 = To_xyz * Tn_o;

% new lattice parameters
obj.lattice.lat_const = sqrt(sum(basisvector2.^2,1));
% new angles
bnorm = bsxfun(@rdivide,basisvector2,obj.lattice.lat_const);
obj.lattice.angle = acos(sum([bnorm(:,2).*bnorm(:,3) bnorm(:,1).*bnorm(:,3) bnorm(:,1).*bnorm(:,2)],1));

% coordinates of the corners of the new coordinate system
pp = [zeros(3,1) Tn_o Tn_o(:,1)+Tn_o(:,3) Tn_o(:,2)+Tn_o(:,3) Tn_o(:,1)+Tn_o(:,2) sum(Tn_o,2)];

% number of cells needed for the extension
nExt  = ceil(max(pp,[],2) - min(pp,[],2))'+2;

% generated atoms
atomList   = obj.atom;
% original number of atoms in the unit cell
nAtom0 = numel(atomList.idx);
atomList.S = obj.unit_cell.S(atomList.idx);
atomList = sw_extendlattice(nExt,atomList);
rExt   = bsxfun(@plus,bsxfun(@times,atomList.RRext,nExt'),(min(pp,[],2)-1));
idxExt = atomList.idxext;


rExt = bsxfun(@plus,rExt,bshift);
% atomic positions in the new unit cell
rNew = inv(Tn_o)*rExt; %#ok<MINV>

epsilon = 10*eps;
% cut atoms outside of the unit cell
%idxCut = any((rNew<0) | (rNew>=(1-eps)),1);
idxCut = any((rNew<-epsilon) | (rNew>(1-epsilon)),1);
rNew(:,idxCut) = [];
idxExt(idxCut) = [];
atomList.Sext(idxCut) = [];

% atoms are close to the face or origin --> put them exactly
rNew(rNew<epsilon) = 0;

% new lattice simmetry is no-symmetry
obj.lattice.sym     = zeros(3,4,0);
obj.lattice.label   = 'P 0';
% no symmetry operations
obj.sym = false;

% new unit cell defined
obj.unit_cell.r     = rNew;
obj.unit_cell.S     = atomList.Sext;

obj.unit_cell.ff  = obj.unit_cell.ff(:,:,idxExt);

fNames = fieldnames(obj.unit_cell);
fNames = setdiff(fNames,{'r' 'S' 'ff'});

for ii = 1:numel(fNames)
    obj.unit_cell.(fNames{ii}) = obj.unit_cell.(fNames{ii})(:,idxExt);
end

% reset the magnetic structure and the bonds
obj.initfield({'coupling' 'mag_str'});

% correct for formula units
obj.unit.nformula = obj.unit.nformula*numel(obj.unit_cell.S)/nAtom0;

% transformation from the original reciprocal lattice into the new
% reciprocal lattice
if nargout>0
    varargout{1} = Tn_o';
end

if param.keepq
    % keep the new coordinate system
    obj.unit.qmat = Tn_o';
end

end