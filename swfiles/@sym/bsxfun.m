function C = bsxfun(op,A,B)
% C = BSXFUN(@op,A,B) extends the bsxfun to sym class variables, for any
% other input type it calls the standrd bsxfun function.
%
% See also sym, syms.
%

if ~isa(B,'sym')
    B = sym(B);
end

nA = ndims(A);
nB = ndims(B);

nD = max(nA,nB);

sA = ones(1,nD);
sB = ones(1,nD);

sA(1:nA) = size(A);
sB(1:nB) = size(B);

if ~all(sA) || ~all(sB)
    C = A;
else
    % repeat dimensions of A where B is non-one, A is one
    rA = (sA==1) & (sB~=1);
    rB = (sB==1) & (sA~=1);
    
    sA(~rB) = 1;
    sB(~rA) = 1;
    
    C = op(repmat(A,sB),repmat(B,sA));
end

end