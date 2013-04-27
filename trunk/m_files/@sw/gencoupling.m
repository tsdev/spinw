function obj = gencoupling(obj, varargin)
% generates the COUPLING property of sw object
%
% obj = GENCOUPLING(obj, 'option1', value, ...)
%
% It calculates the shortest distances between magnetic ions and sort them
% according to increasing distance. Two couplings are regarded to be
% equivalent if the interion distances are equal within a defined
% tolerance. This function should be called after the lattice and atomic
% positions are defined.
%
% Options:
%
% nUnitCell     Edge length of the parallelepiped withing the
%               algorithm searches for neares neighbours in lattice
%               units. Default is 3.
% maxDistance   Maximum inter-ion distance that will be stored in the
%               obj.coupling property in units of Angstrom. Default
%               is 6.
% tolDistance   Tolerance of distance, within two couplings are regarded
%               equivalent, default is 1e-5.
%
% See also SW.
%

inpForm.fname  = {'nUnitCell' 'maxDistance' 'tolDistance'};
inpForm.defval = {3           6             1e-5         };
inpForm.size   = {[1 1]       [1 1]         [1 1]        };

param = sw_readparam(inpForm, varargin{:});

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
    
    coupling.dist    = zeros(1,nDist);
    coupling.dl      = zeros(3,nDist);
    coupling.atom1   = zeros(1,nDist);
    coupling.atom2   = zeros(1,nDist);
    coupling.mat_idx = zeros(3,nDist);
    coupling.idx     = zeros(1,nDist);
    
    
    % Atomic position [Angstrom] of the magnetic atoms.
    atcoords       = obj.basisvector*atPos;
    basis_atcoords = obj.basisvector*mAtom.r;
    
    % Inside the original unit cell magn. atoms couple to atoms outside unit cell.
    index=1;
    for ii=1:nMagAtom
        for jj=1:nHalfCube
            coupling.dist(index)  = norm(basis_atcoords(:,ii)-atcoords(:,jj));
            coupling.dl(:,index)  = atTr(:,jj);
            coupling.atom1(index) = ii;
            coupling.atom2(index) = atIndex(jj);
            index = index+1;
        end
    end
    
    % Couplings inside origin unit cell.
    for ii=1:(nMagAtom-1)
        for jj=(ii+1):nMagAtom
            coupling.dist(index)  = norm(basis_atcoords(:,ii)-basis_atcoords(:,jj));
            coupling.dl(:,index)  = [0 0 0];
            coupling.atom1(index) = ii;
            coupling.atom2(index) = jj;
            index = index+1;
        end
    end
    
    [~, newindex]  = sort(coupling.dist);
    coupling.dist  = coupling.dist(newindex);
    coupling.dl    = coupling.dl(:,newindex);
    coupling.atom1 = coupling.atom1(newindex);
    coupling.atom2 = coupling.atom2(newindex);
    
    % Clear couplings larger than maxDistance [Angstrom].
    index = coupling.dist<param.maxDistance;
    coupling.dist    = coupling.dist(index);
    coupling.dl      = int32(coupling.dl(:,index));
    coupling.atom1   = int32(coupling.atom1(index));
    coupling.atom2   = int32(coupling.atom2(index));
    coupling.idx     = int32(coupling.idx(index));
    coupling.mat_idx = int32(coupling.mat_idx(:,index));
    
    % New number of couplings.
    nDist = size(coupling.dist,2);
    
    % Finds the equivalent distances.
    runVal = coupling.dist(1);
    
    index = 1;
    for ii=1:nDist
        if abs(coupling.dist(ii)-runVal)>param.tolDistance
            index = index+1;
            runVal = coupling.dist(ii);
        end
        coupling.idx(ii) = index;
    end
    
    aniso = int32(zeros(1,nMagAtom));
    
else
    % If there is no magnetic atom the coupling and anisotropy are empty.
    coupling.dist    = [];
    coupling.dl      = int32(zeros(3,0));
    coupling.atom1   = int32(zeros(1,0));
    coupling.atom2   = int32(zeros(1,0));
    coupling.mat_idx = int32(zeros(3,0));
    coupling.idx     = int32(zeros(1,0));
    
    aniso            = int32(zeros(1,0));
    
end

coupling             = rmfield(coupling,'dist');
obj.coupling         = coupling;
obj.single_ion.aniso = aniso;

validate(obj);

end