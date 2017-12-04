function L = sw_bonddim(C, nAtom)
% find dimensionality of a periodic bond network
% 
% ### Syntax
% 
% `l = sw_bonddim(c, {natom})`
% 
% ### Description
% 
% `l = sw_bonddim(c, {natom})` splits the given periodic bond network into
% disjunct subsystems and determines the dimensionality of each subsystem.
% 
% ### Examples
%
% Check the bond dimensionality of the triangular lattice:
%
% ```
% >>tri = sw_model('triAF')>>
% >>sw_bonddim(tri.intmatrix.all)>>
% ```
%
% ### Input Arguments
% 
% `C`
% : Bond list in a matrix with dimensions of $[5\times n_{bond}]$, where the meaning of
%   the rows are:
%   * `#1:#3`   Lattice translations between the coupled atoms in
%               lattice units (always integer).
%   * `#4`      Index of the bond starting atom.
%   * `#5`      Index of the bond end atom.
%
%   For example for a chain along b-axis on a Bravais lattice:
%       `C = [1;1;0;1;0]`
% 
% `nAtom`
% : Number of atoms in the unit cell. If not given, the maximum
%   atom index from the bond list is taken.
% 
% ### Output Arguments
% 
% `L`
% : Struct with the number of elements equal to the number of
%           subsystems, it has the following fields:
%   * `D`       Dimensionality of the subsystem $(0\leq D\leq 3)$.
%   * `base`    Basis vectors spanning the subsystem stored in a
%                       $[3\times D]$ matrix where each column denotes a basis
%                       vector.
%   * `site`    List of sites that belong to the subsystem.
%

% reshuffle C to match the input with the output of spinw.intmatrix
% in the code
%       C(1,:)      is atom1
%       C(2,:)      is atom2
%       C(3:5,:)    is the transplation vector
C = C([4:5 1:3],:);

% input check
if size(C,1)~=5
    error('sw_bonddim:WrongInput','The given bond matrix has wrong dimensions!')
end
% all values have to be integer
if any(floor(C(:))-ceil(C(:)))
    error('sw_bonddim:WrongInput','The given bond matrix must contain only integers!')
end

% check that all index up to a maximum is present
atom = [C(1,:) C(2,:)];

if nargin<2
    nAtom = max(atom);
end

% add missing indices
missing = find(~ismember(1:nAtom,atom));
if ~isempty(missing)
    C = [C [missing;missing;zeros(3,numel(missing))]];
end

% sort bonds in increasing order
idx = C(2,:)<C(1,:);
C(:,idx) = flip(C(:,idx));

% store subsystems in struct
L = struct('D',cell(1,0),'base',cell(1,0),'site',cell(1,0));
iSub = 1;

% find other loops going through every subsystem
while ~isempty(C)
    % start with lowest index atom as origin
    atom0 = min(C(1,:));
    
    % starting atom
    L0 = [atom0;zeros(3,1)];
    
    while any(ismember(C(1,:),L0(1,:)) | ismember(C(2,:),L0(1,:)))
        % find self loops
        idx0 = ismember(C(1,:),L0(1,:)) & C(1,:)==C(2,:);
        % find next neighbors
        idx1 = find(ismember(C(1,:),L0(1,:)) & ~idx0);
        idx2 = find(ismember(C(2,:),L0(1,:)) & ~idx0);
        C(:,idx2) = flip(C(:,idx2));
        idx = [find(idx0) unique([idx1 idx2])];
        
        % add new atoms to the L list
        for ii = 1:numel(idx)
            dL = bsxfun(@plus,L0(2:4,L0(1,:)==C(1,idx(ii))),C(3:5,idx(ii)));
            L0  = [L0 [repmat(C(2,idx(ii)),1,size(dL,2));dL]]; %#ok<AGROW>
        end
        % remove neighbors from bond list
        C(:,idx) = [];
    end
    
    % find all combination of translation vectors within subsystem
    idx  = unique(L0(1,:));
    % basis vectors
    V = zeros(3,0);
    ii = 1;
    while ii <= numel(idx) && size(V,2)<3
        lSel = L0(2:4,L0(1,:)==idx(ii));
        % normalize vectors
        %lSel = bsxfun(@rdivide,lSel,sum(lSel.^2,1));
        %lSel(isnan(lSel)) = 0;
        % keep unique ones only
        %lSel = sw_uniquetol(lSel,1e-6);
        lSel = unique(lSel','rows')';
        V = orth([V reshape(bsxfun(@minus,permute(lSel,[2 3 1]),permute(lSel,[3 2 1])),[],3)']);
        ii = ii + 1;
    end
    
    % store the subsystem basis vectors
    L(iSub).base = V;
    L(iSub).D    = size(V,2);
    if L(iSub).D == 3
        % pretty 3D
        L(iSub).base = eye(3);
    end
    L(iSub).site = unique([atom0 idx]);
    iSub = iSub + 1;
end

% sort according to decreasing dimensionality
[~,idx] = sort([L.D],'descend');
L = L(idx);

end

function C = flip(C)
% flip bond atom1 --> atom2

C = [C(2,:);C(1,:);-C(3:5,:)];

end