function [V,D] = eigorth(M, tol)
% orthogonal eigenvectors of defective eigenvalues
% [V, D] = EIGORTH(M, tol)
%
% Input:
%
% M     Square matrix.
% tol   Tolerance, within two eigenvalue are regarded equal:
%           abs(real(diff(D))) < max([real(D) 1e-5])*tol
%       Optional.
%
% Output:
%
% V     Matrix, every column contains an eigenvector of M.
% D     Vector, that contains the eigenvalues of M.
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

[V, D] = eig(M);

D = diag(D)';

% sort eigenvalues-eigenvectors according to the real part of D
[~, idx] = sort(real(D));
D = D(idx);
V = V(:,idx);
% selects degenerate eigenvalues
degIdx = abs(real(diff(D))) < max([real(D) 1e-5])*tol;
degIdx = [0 degIdx];
% degIdx contains islands with increasing integers
degVal = cumsum(diff([0 degIdx])==1);
degVal = [degVal(2:end) degVal(end)];
degIdx = (degIdx | [degIdx(2:end) 0]).*degVal;

try
    % orthogonalise them
    for ii = 1:max(degIdx)
        V(:,degIdx==ii) = orth(V(:,degIdx==ii));
    end
catch %#ok<CTCH>
    % if there are not enough orthogonal vector, skip the orthogonalisation
    warning('eigorth:Error','Eigenvectors of defective eigenvalues cannot be orthogonalised!');
    [V, D] = eig(M);
    
    D = diag(D)';
end
end