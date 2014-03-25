function varargout = newcell(obj,bvect, bshift)
% changes lattice vectors while keeping atoms
%
% {T} = NEWCELL(obj, bvect, {bshift}) defines new unit cell using 3 vectors
% contained in bvect cell: {v1, v2, v3}. The three vectors in lattice units
% define a parallelepiped. This will be the new unit cell. The atoms from
% the original unit cell will fill the new unit cell. bshift vector defines
% a shift of the position of the unit cell.
% There might be problem, if there are atoms on the faces and corners of
% the new lattice (due to numerical error). In this case use small bshift
% to move the atoms.
% 
% Output:
%
% T     is a transformation matrix that converts Q points (in reciprocal
%       lattice units) from the old reciprocal lattice to the new
%       reciprocal lattice as follows:
%           Qrlu_new = T * Qrlu_old,
%       where the dimensions of the Q vectors are [1 3].
%
% See also SW.GENLATTICE.
%

if nargin <= 1
    help sw.newcell;
    return;
end

if ~iscell(bvect) || numel(bvect)~=3
    error('sw:newcell:WrongInput','Input has to be cell type with 3 vectors inside!');
end

% shift
if nargin == 2
    bshift = [0;0;0];
else
    bshift = bshift(:);
end

% here 3 coordinate systems are used:
% - xyz real space, Cartesian
% - original unit cell a,b and c vectors
% - new unit cell a',b' and c' vectors


% transformation matrix from the new lattice units into the original
% lattice
% v_orig = Tn_o*v_new
Tn_o = [bvect{1}(:) bvect{2}(:) bvect{3}(:)];

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
%obj.mag_str.N_ext = int32(nExt);

% generated atoms
atomList = obj.atom;
atomList.S = obj.unit_cell.S(atomList.idx);
atomList = sw_extendlattice(nExt,atomList);
rExt   = bsxfun(@plus,bsxfun(@times,atomList.RRext,nExt'),(min(pp,[],2)-1));
idxExt = atomList.idxext;
Sext   = atomList.Sext;

rExt = bsxfun(@plus,rExt,bshift);
% atomic positions in the new unit cell
rNew = inv(Tn_o)*rExt; %#ok<MINV>

% cut atoms outside of the unit cell
idxCut = any((rNew<0) | (rNew>=1),1);
rNew(:,idxCut) = [];
idxExt(idxCut) = [];
Sext(idxCut)   = [];

% new lattice simmetry is no-symmetry
obj.lattice.sym     = int32(0);
% no symmetry operations
obj.sym = false;

% new unit cell defined
obj.unit_cell.r     = rNew;
obj.unit_cell.S     = Sext;
obj.unit_cell.label = obj.unit_cell.label(idxExt);
obj.unit_cell.color = obj.unit_cell.color(:,idxExt);
obj.mag_str.N_ext   = int32([1 1 1]);

% transformation from the original reciprocal lattice into the new
% reciprocal lattice
if nargout>0
    varargout{1} = inv(Tn_o);
end

end