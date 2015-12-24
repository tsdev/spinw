function [V,D] = eigorth(M, tol, sortMode, useMex)
% orthogonal eigenvectors of defective eigenvalues
% [V, D, {Idx}] = EIGORTH(M, tol, sort)
%
% Input:
%
% M         Stack of square matrices, dimensions are [nMat nMat nStack].
% tol       Tolerance, within two eigenvalue are regarded equal:
%               abs(real(diff(D))) < max([real(D) 1e-5])*tol
%           Default is 1e-5.
% sortMode  If true, mode sorting is performed, and idx stores the sorted
%           permutation of the modes and eigenvectors. Default is false.
%
% Output:
%
% V         Matrix stack, every column contains an eigenvector of M in
%           every stack, dimensions are the same as M.
% D         Stack of vectors, that contains the eigenvalues of M,
%           dimensions are [nMat nStack].
% idx       Permutation indices of the eigenvalues and eigenvectors for
%           every stack, dimensions are [nMat nStack].
%
% See also eig.
%

if nargin == 0
    help eigorth;
    return;
end

if nargin < 2
    tol = 1e-5;
end

if nargin < 3
    sortMode = false;
end

if nargin < 4
    useMex = false;
end

nMat   = size(M,1);
nStack = size(M,3);

% Use OpenMP parallelised mex file if it exists
if nStack>1 && useMex && exist('eig_thr')==3
    if sortMode
        [V, D] = eigenshuffle(M,useMex,'orth');
    else
        [V, D] = eig_thr(M,'orth');
    end
    return;
end
    
if sortMode
    [V, D] = eigenshuffle(M);
    
else
    V = zeros(nMat,nMat,nStack);
    D = zeros(nMat,nStack);
    
    for ii = 1:nStack
        [V(:,:,ii), Dtemp] = eig(M(:,:,ii));
        D(:,ii) = diag(Dtemp);
    end
end

for jj = 1:nStack
    % sort eigenvalues-eigenvectors according to the real part of D
    [~, idxD] = sort(real(D(:,jj)));
    Di = D(idxD,jj)';
    Vi = V(:,idxD,jj);
    % selects degenerate eigenvalues
    degIdx = [0 abs(real(diff(Di))) < max([real(Di) 1e-5])*tol];
    % degIdx contains islands with increasing integers
    degVal = cumsum(diff([0 degIdx])==1);
    degVal = [degVal(2:end) degVal(end)];
    degIdx = (degIdx | [degIdx(2:end) 0]).*degVal;
    
    try
        % orthogonalise them
        for ii = 1:max(degIdx)
            Vi(:,degIdx==ii) = orth(Vi(:,degIdx==ii));
        end
        % permute back Vi
        Vi(:,idxD) = Vi;
        
        V(:,:,jj) = Vi;
    catch %#ok<CTCH>
        % if there are not enough orthogonal vector, skip the orthogonalisation
        warning('eigorth:Error','Eigenvectors of defective eigenvalues cannot be orthogonalised!');
    end
end

end
