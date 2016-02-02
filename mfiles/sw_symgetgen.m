function [symOpG, symTrG, isGen] = sw_symgetgen(symOp, symTr)
% creates the generators from a list of symmetry operators
%
% [symOp, symTr, isGen] = SW_SYMGETGEN(symOp, symTr)
%
% Input:
%
% symOp     Matrix, containing the rotational part of the symmetry
%           operators, dimensions are [3 3 nSym].
% symTr     Matrix, whose columns contain the translation of every symmetry
%           operator, dimensions are [3 nSym].
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

uIdx = 1:size(symOp,3);
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
    

while ~isempty(uIdx)
    symOpG = cat(3,symOpG,symOp(:,:,uIdx(1)));
    symTrG = [symTrG symTr(:,uIdx(1))]; %#ok<AGROW>
    
    % generate all symmetry elements from the given generators
    [symOp1, symTr1] = sw_gencoord({symOpG symTrG});
    % number of already generated symmetry operators
    nReady = size(symOp1,3);
    % find all non-generated operators
    symMat = [[reshape(symOp1,9,[]) reshape(symOp,9,[])];[symTr1 symTr]]';
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
    symTrG = zeros(1,3);
    isGen  = 1;
end

end