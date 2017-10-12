function L = sw_finddim3(C)
% find dimensionality of a periodic bond network
%
% L = SW_BONDDIM(C)
%
% The function splits the given periodic bond network into disjunct
% subsystems and determines the dimensionality of each subsystem.
%
% Input:
%
% C         Bond list in a matrix with dimensions [5 nBond]. The meaning of
%           the rows:
%               #1      Index of the bond starting atom.
%               #2      Index of the bond end atom.
%               #3:#5   Lattice translations between the coupled atoms in
%                       lattice units (always integer).
%           For example for a chain along b-axis on a Bravais lattice:
%               C = [1;1;0;1;0]
%
% Output:
%
% L         Struct with the number of elements equal to the number of
%           subsystems, it has the following fields:
%               D       Dimensionality of the subsystem (0<=D<=3).
%               base    Basis vectors spanning the subsystem stored in a
%                       [3 D] matrix where each clumn denotes a basis
%                       vector.
%               site    List of sites that belong to the subsystem.
%


% input check
if size(C,1)~=5
    error('sw_finddim:WrongInput','The given bond matrix has wrong dimensions!')
end
% all values have to be integer
if any(floor(C(:))-ceil(C(:)))
    error('sw_finddim:WrongInput','The given bond matrix must contain only integers!')
end
% check that all index up to a maximum is present
atom = [C(1,:) C(2,:)];
% add missing indices
missing = find(~ismember(1:max(atom),atom));
C = [C [missing;missing;zeros(3,numel(missing))]];

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
    allL = zeros(3,0);
    for ii = 1:numel(idx)
        lSel = L0(2:4,L0(1,:)==idx(ii));
        allL = [allL reshape(bsxfun(@minus,permute(lSel,[2 3 1]),permute(lSel,[3 2 1])),[],3)']; %#ok<AGROW>
    end
    % store the subsystem basis vectors
    V = orth(allL);
    L(iSub).base = V;
    L(iSub).D    = size(V,2);
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