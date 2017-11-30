function [unique, firstIdx] = sw_uniquetol(M,tol)
% returns the unique column vectors within tolerance
% 
% ### Syntax
% 
% `[Mu, firstIdx] = sw_uniquetol(M,tol)`
% 
% ### Description
% 
% `[Mu, firstIdx] = sw_uniquetol(m,tol)` returns unique column vectors
% within the given `tol` tolerance. Two column vectors are considered
% unequal, if the distance between them is larger than the tolerance
% ($\delta$):
%
% $\sqrt{\sum_i (V_i-U_i)^2} < \delta$
% 
% ### Input Arguments
% 
% `M`
% : Matrix that contains column vectors.
% 
% `tol`
% : Distance tolerance, default value is $10^{-5}$.
% 
% ### Output Arguments
% 
% `Mu`
% : Matrix that contains the unique column vectors.
%
% `firstIdx`
% : Indices pointing to the first occurence of the unique element.
%
% This function is similar to the Matlab built-in
% `unique(M,'rows','first')`, but with controllable tolerance.
%
% ### See Also
%
% [unique](https://ch.mathworks.com/help/matlab/ref/unique.html)
%

if nargin == 0
    swhelp sw_uniquetol
    return
end

if nargin < 2
    tol = 1e-5;
end

unique = M*0;
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
        idxSame = sum(bsxfun(@minus,M,unique(:,idx)).^2,1) < tol2;
        M(:,idxSame)  = [];
        idxM(idxSame) = [];
        
        idx = idx + 1;
    end
    firstIdx = firstIdx(1:(idx-1));
end
unique = unique(:,1:(idx-1));

end