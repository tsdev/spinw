function [unique, firstIdx] = sw_uniquetol(M,tol)
% returns the unique column vectors within tolerance
%
% [unique, firstIdx] = SW_UNIQUETOL(M,tol)
%
% Two column vectors are considered unequal, if the distance between them
% is larger than the tolerance.
%
% Input:
%
% M         Matrix that contains column vectors.
% tol       Distance tolerance, default is 1e-5.
%
% Output:
%
% unique    Matrix that contains the unique column vectors.
% firstIdx  Indices pointing to the first occurence of the unique element.
%
% This function is similar to the Matlab built-in unique(M,'rows','first'),
% but with arbitrary tolerance.
%

if nargin == 0
    help sw_uniquetol
    return
end

if nargin < 2
    tol = 1e-5;
end

unique = zeros(size(M));
tol2 = tol(1)^2;

if nargout < 2
    idx = 1;
    while ~isempty(M)
        unique(:,idx) = M(:,1);
        idxSame = sum(bsxfunsym(@minus,M,unique(:,idx)).^2,1) < tol2;
        M(:,idxSame) = [];
        idx = idx + 1;
    end
    
else
    idx  = 1;
    % storing the indices in M
    idxM = 1:size(M,2);
    firstIdx = zeros(1,size(M,2));
    
    while ~isempty(M)
        unique(:,idx) = M(:,1);
        firstIdx(idx) = idxM(1);
        idxSame = sum(bsxfunsym(@minus,M,unique(:,idx)).^2,1) < tol2;
        M(:,idxSame)  = [];
        idxM(idxSame) = [];
        
        idx = idx + 1;
    end
    firstIdx = firstIdx(1:(idx-1));
end
unique = unique(:,1:(idx-1));

end