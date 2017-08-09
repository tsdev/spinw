function Lsub = sw_finddim(C)
% find dimensionality of a periodic bond network

% inpu check
if size(C,1)~=5
    error('sw_finddim:WrongInput','The given bond matrix has wrong dimensions!')
end
% check that all index up to a maximum is present
idxList = unique([C(1,:) C(2,:)]);
if any(floor(C(:))-ceil(C(:)))
    error('sw_finddim:WrongInput','The given bond matrix must contain only integers!')
end
if min(idxList)~=1 || max(idxList)~=numel(idxList)
    error('sw_finddim:WrongInput','The given bond matrix does not contain consecutive atom indices!')
end


% find self loops
idx = find(C(1,:)==C(2,:));
% L keeps atom index and translation from origin
L0 = C([1 3:5],idx);
C(:,idx) = [];

% sort bonds in increasing order
idx = C(2,:)<C(1,:);
C(:,idx) = flip(C(:,idx));

% store subsystems in struct
Lsub = struct('base',cell(1,0),'site',cell(1,0));

% find other loops going through every subsystem
while ~isempty(C)
    % start with lowest index atom as origin
    atom0 = min(C(1,:));
    idx   = find(C(1,:) == atom0);
    % store the first neighbors
    L = C(2:end,idx);
    C(:,idx) = [];
    
    % find 2nd neighbors
    idx1 = find(ismember(C(1,:),L(1,:)));
    idx2 = find(ismember(C(2,:),L(1,:)));
    C(:,idx2) = flip(C(:,idx2));
    idx = [idx1 idx2];
    
    while ~isempty(idx)
        atom1 = C(1,idx);
        atom2 = C(2,idx);
        % add new atoms to the L list
        for ii = 1:numel(idx)
            dL = bsxfun(@plus,L(2:4,L(1,:)==atom1(ii)),C(3:5,idx(ii)));
            L  = [L [repmat(atom2(ii),1,size(dL,2));dL]]; %#ok<AGROW>
        end
        % remove neighbors
        C(:,idx) = [];
        % n+1 neighbors
        % find neighbors
        idx1 = find(ismember(C(1,:),L(1,:)));
        idx2 = find(ismember(C(2,:),L(1,:)));
        C(:,idx2) = flip(C(:,idx2));
        idx = [idx1 idx2];
    end
    
    % find all combination of translation vectors within subsystem
    idx  = unique(L(1,:));
    allL = L0(2:4,ismember(L0(1,:),idx));
    for ii = 1:numel(idx)
        lSel = L(2:4,L(1,:)==idx(ii));
        allL = [allL reshape(bsxfun(@minus,permute(lSel,[2 3 1]),permute(lSel,[3 2 1])),[],3)']; %#ok<AGROW>
    end
    % store the subsystem basis vectors
    V = orth(allL);
    Lsub(end+1).base = {V}; %#ok<AGROW>
    Lsub(end).dim    = size(V,2);
    Lsub(end).site   = [atom0 idx];
end

% add any remaining self loops

end

function C = flip(C)
% flip bond atom1 --> atom2

C = [C(2,:);C(1,:);-C(3:5,:)];

end