function sumA = symsum2(A,dim)

nA = ndims(A);
dimA = 1:nA;
dimA(1)  = dim;
dimA(dim) = 1;

A = permute(A,dimA);
sA = size(A);

sA(1) = 1;


sumA = A(1,:);

if dA > 1
    for ii = 2:dA
        sumA = sumA + A(ii,:);
    end
end

sumA = reshape(sumA,sA);

sumA = permute(sumA,dimA);

end