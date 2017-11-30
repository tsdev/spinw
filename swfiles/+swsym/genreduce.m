function [symOpG, isGen] = genreduce(symOp)
% reduces symmetry operators to the generators
% 
% ### Syntax
% 
% `[symOpG, isGen] = swsym.genreduce(symOp)`
% 
% ### Description
% 
% `[symOpG, isGen] = swsym.genreduce(symOp)` takes the list of symmetry
% operators in `symOp` and determines a minimum subset of operators that
% can generate the given list.
% 
% ### Input Arguments
% 
% `symOp`
% : Matrix that contains both the rotation and translation matrices
%   having dimensions of $[3\times 4\times n_{sym}]$, where the
%   `symMat(:,4,:)` stores the translation vectors, while the
%   `symMat(:,1:3,:)` stores the $3\times 3$ rotation matrices.
% 
% ### Output Arguments
% 
% `symOpG`
% : A set of operators, that can generate all the operators of the input.
%
% `isGen`
% : Vector, that gives whether a given input operator is part of
%   the generators, dimensions are $[1\times n_{sym}]$.
% 
% ### See Also
% 
% [swsym.add] \| [swsym.generator] \| [swsym.operator]
%

if nargin == 0
    swhelp swsym.genreduce
    return
end

symOpG = zeros(3,4,0);

% original list of operators
symOp0 = symOp;

% find pure translation operator and put to the first place
pureT = find(~sumn(abs(bsxfun(@minus,symOp(:,1:3,:),eye(3))),[1 2]) & sum(symOp(:,4,:)~=0,1));
symOp = symOp(:,:,[pureT' 1:end]);

while ~isempty(symOp)
    symOpG = cat(3,symOpG,symOp(:,:,1));
    % generate all symmetry elements from the given generators
    symOp1 = swsym.operator(symOpG);
    symOp  = reshape(setdiff(reshape(symOp,12,[])',reshape(symOp1,12,[])','rows','stable')',3,4,[]);
end

% remove any additional operator that is not necessary
idx = 1;
while idx<=size(symOpG,3)
    symSel = symOpG(:,:,idx);
    symOpG(:,:,idx) = [];
    if ~isempty(setdiff(reshape(symSel,12,[])',reshape(swsym.operator(symOpG),12,[])','rows'))
        symOpG = cat(3,symSel,symOpG);
        idx = idx+1;
    end
end

% determine isGen
isGen = ismember(reshape(symOp0,12,[])',reshape(symOpG,12,[])','rows');

% P1
if isempty(symOpG)
    symOpG = [eye(3) zeros(3,1)];
    isGen  = 1;
end

end