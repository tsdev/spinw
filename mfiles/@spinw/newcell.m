function varargout = newcell(obj,bvect, bshift)
% changes lattice vectors while keeping atoms
%
% {T} = NEWCELL(obj, bvect, {bshift}) 
%
% The function defines new unit cell using the 3 vectors contained in
% bvect. The three vectors in lattice units define a parallelepiped. This
% will be the new unit cell. The atoms from the original unit cell will
% fill the new unit cell.
%
% Input:
%
% obj       spinw class object.
% bvect     Defines the new lattice vectors in the original lattice
%           coordinate system. Cell with the following elements 
%           {v1 v2 v3}.
% bshift    Vector defines a shift of the position of the unit cell.
%           Optional.
%
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
% Example:
%
% tri = sw;
% tri.genlattice('lat_const',[3 3 5],'angled',[90 90 120])
% tri.addatom('r',[0 0 0])
% tri.newcell({[1 0 0] [1 2 0] [0 0 1]})
% plot(tri)
%
% The example show how to convert a triangular lattice into orthogonal
% lattice vectors and plots the new unit cell.
%
% See also SPINW.GENLATTICE, SPINW.GENCOUPLING, SPINW.NOSYM.
%

if nargin <= 1
    help spinw.newcell
    return
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
atomList   = obj.atom;
atomList.S = obj.unit_cell.S(atomList.idx);
atomList = sw_extendlattice(nExt,atomList);
rExt   = bsxfun(@plus,bsxfun(@times,atomList.RRext,nExt'),(min(pp,[],2)-1));
idxExt = atomList.idxext;



rExt = bsxfun(@plus,rExt,bshift);
% atomic positions in the new unit cell
rNew = inv(Tn_o)*rExt; %#ok<MINV>

% cut atoms outside of the unit cell
%idxCut = any((rNew<0) | (rNew>=(1-eps)),1);
idxCut = any((rNew<0) | (rNew>=1),1);
rNew(:,idxCut) = [];
idxExt(idxCut) = [];
atomList.Sext(idxCut)   = [];

% find equivalent positions for dubious border atoms
idxB = find(any(rNew>=(1-eps),1));
% find the axis which is problematic
idxA = rNew(:,idxB)>=(1-eps);
% loop over and check whether the atom exists already
idxCut = [];
for ii = 1:numel(idxB)
    rTemp = rNew(:,idxB(ii));
    rTemp(idxA(:,ii)) = 0;
    % atoms indices in rNew that are equivalent to the bad atom
    if any(sum(bsxfun(@minus,rNew,rTemp).^2,1)<eps)
        % remove the bad behaving atom
        idxCut = [idxCut idxB(ii)]; %#ok<AGROW>
    end
end

% remove the necessary bad behaving atoms
rNew(:,idxCut) = [];
idxExt(idxCut) = [];
atomList.Sext(idxCut)   = [];

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

% obj.unit_cell.label = obj.unit_cell.label(:,idxExt);
% obj.unit_cell.color = obj.unit_cell.color(:,idxExt);
% obj.unit_cell.ox    = obj.unit_cell.ox(:,idxExt);
% obj.unit_cell.occ   = obj.unit_cell.occ(:,idxExt);
% obj.unit_cell.b     = obj.unit_cell.b(:,idxExt);
% obj.unit_cell.A     = obj.unit_cell.A(:,idxExt);
% obj.unit_cell.Z     = obj.unit_cell.Z(:,idxExt);
% obj.unit_cell.biso  = obj.unit_cell.biso(:,idxExt);


obj.mag_str.N_ext   = int32([1 1 1]);


% transformation from the original reciprocal lattice into the new
% reciprocal lattice
if nargout>0
    varargout{1} = Tn_o;
end

end