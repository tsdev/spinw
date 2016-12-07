function [symOpG, symTrG, isGen] = sw_symgetgen(symOp, symTr)
% creates the generators from a list of symmetry operators
%
% [symOp, symTr, isGen] = SW_SYMGETGEN(symOp, symTr)
% [symOp, symTr, isGen] = SW_SYMGETGEN(symMat)
%
% Input:
%
% symOp     Matrix, containing the rotational part of the symmetry
%           operators, dimensions are [3 3 nSym].
% symTr     Matrix, whose columns contain the translation of every symmetry
%           operator, dimensions are [3 nSym].
% symMat    Matrix that contains both the rotation and translation matrices
%           having dimensions of [3 4 nSym], where the symMat(:,4,:) stores
%           the translation vectors, while the symMat(:,1:3,:) stores the
%           rotation operators.
%
% Output:
%
% symOp, symTr
%           A set of operators, that can generate all the operators of the
%           input.
% isGen     Vector, that gives whether a given input operator is part of
%           the generators, dimensions are [1 nSym].
%

if nargin == 0
    help sw_symgetgen
    return
end

if nargin == 1
    % use the general input
    symTr = permute(symOp(:,4,:),[1 3 2]);
    symOp = symOp(:,1:3,:);
end

uIdx   = 1:size(symOp,3);
symOpG = zeros(3,3,0);
symTrG = zeros(3,0);

% For Matlab earlier than R2012a, the unique function doesn't need the
% option R2012a and older versions doesn't recognise the 'R2012a' flag,
% they throw an error.
if verLessThan('matlab','7.14')
    uOption = {'rows'};
else
    uOption = {'rows','R2012a'};
end

% run twice
symOp0 = symOp;
symTr0 = symTr;

% eliminate operators twice
for ii = 1:2
    while ~isempty(uIdx)
        symOpG = cat(3,symOpG,symOp0(:,:,uIdx(1)));
        symTrG = [symTrG symTr0(:,uIdx(1))]; %#ok<AGROW>
        
        % generate all symmetry elements from the given generators
        [symOp1, symTr1] = sw_gencoord({symOpG symTrG});
        % number of already generated symmetry operators
        nReady = size(symOp1,3);
        % find all non-generated operators
        symMat = [[reshape(symOp1,9,[]) reshape(symOp0,9,[])];[symTr1 symTr0]]';
        [~, uIdx] = unique(symMat,uOption{:});
        uIdx = sort(uIdx)';
        uIdx = uIdx(nReady+1:end)-nReady;
    end
    
    % remove the unity operator
    symMat = [[reshape(eye(3),9,[]) reshape(symOpG,9,[])];[[0;0;0] symTrG]]';
    [~, uIdx] = unique(symMat,uOption{:});
    uIdx = sort(uIdx)';
    uIdx = uIdx(2:end)-1;
    
    symOpG = symOpG(:,:,uIdx);
    symTrG = symTrG(:,uIdx);
    
    if ii == 1
        % reshuffle the operators
        symOp0 = symOpG(:,:,end:-1:1);
        symTr0 = symTrG(:,end:-1:1);
        uIdx   = 1:size(symOp0,3);
        symOpG = zeros(3,3,0);
        symTrG = zeros(3,0);
    end
end

% determine isGen
symMat = [[reshape(symOpG,9,[]) reshape(symOp,9,[])];[symTrG symTr]]';
[~, uIdx] = unique(symMat,uOption{:});
uIdx = sort(uIdx)';
nGen = size(symOpG,3);
uIdx = uIdx((nGen+1):end)-nGen;

isGen = true(1,size(symOp,3));
isGen(uIdx) = false;

% P1
if isempty(symOpG)
    symOpG = eye(3);
    symTrG = zeros(3,1);
    isGen  = 1;
end

end