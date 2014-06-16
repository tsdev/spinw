function sumA = sumsym(A,dim)
% sumA = SUMSYM(A,dim) sums up matrices containing symbolic variables in
% arbitrary dimensions, for any other input type it calls the standard sum
% function.
%
% See also sym, syms.
%

if isa(A,'sym')
    nSum = size(A,dim);
    
    nA = ndims(A);
    dimA = 1:nA;
    dimA(1)  = dim;
    dimA(dim) = 1;
    
    A = permute(A,dimA);
    sA = size(A);
    
    sA(1) = 1;
    
    
    sumA = A(1,:);
    
    if  nSum > 1
        for ii = 2:nSum
            sumA = sumA + A(ii,:);
        end
    end
    
    sumA = reshape(sumA,sA);
    
    sumA = permute(sumA,dimA);
    
else
    sumA = sum(A,dim);
end

end