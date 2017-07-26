function D = inv3(M,opt)
% calculates the matrix inverse of 3x3 stacked matrices
%
% D = INV3(M,{opt})
%
% The function is vectorised and efficient for very large number of
% matrices. If opt is 'diag', only the diagonal of the inverse matrices are
% calculated.
%
% Carefull, the code is numerically not stable!!!
%
% Input:
%
% M         Matrix with dimensions of 3 x 3 x N1 x N2 x ...
%
% Output:
%
% D         Matrix with dimensions of 3 x 3 x N1 x N2 x ..., if opt ==
%           'diag', the output is 3 x N1 x N2 x ... matrix.
%

if nargin == 0
    help inv3
    return
end

diagOpt = false;
if nargin > 1 && strcmp(opt,'diag')
    diagOpt = true;
end

sM = [size(M) 1];

if sM(1)~= 3 || sM(2) ~= 3
    error('invdiag:WrongInput','Input matrix has wrong dimensions!')
end

% deal with vectors
M = reshape(M,3,3,[]);

detM = M(1,1,:).*(M(2,2,:).*M(3,3,:)-M(2,3,:).*M(3,2,:))...
    +M(1,2,:).*(M(2,3,:).*M(3,1,:)-M(2,1,:).*M(3,3,:))...
    +M(1,3,:).*(M(2,1,:).*M(3,2,:)-M(2,2,:).*M(3,1,:));

if diagOpt
    D = zeros(3,size(M,3));
    
    D(1,:) = M(2,2,:).*M(3,3,:)-M(2,3,:).*M(3,2,:);
    D(2,:) = M(1,1,:).*M(3,3,:)-M(1,3,:).*M(3,1,:);
    D(3,:) = M(1,1,:).*M(2,2,:)-M(1,2,:).*M(2,1,:);
    
    D = bsxfun(@rdivide,D,permute(detM,[1 3 2]));
    D = reshape(D,[3 sM(3:end)]);
    
else
    D = zeros(3,3,size(M,3));
    
    D(1,1,:) = M(2,2,:).*M(3,3,:)-M(2,3,:).*M(3,2,:);
    D(1,2,:) = M(1,3,:).*M(3,2,:)-M(1,2,:).*M(3,3,:);
    D(1,3,:) = M(1,2,:).*M(2,3,:)-M(1,3,:).*M(2,2,:);
    
    
    D(2,1,:) = M(2,3,:).*M(3,1,:)-M(2,1,:).*M(3,3,:);
    D(2,2,:) = M(1,1,:).*M(3,3,:)-M(1,3,:).*M(3,1,:);
    D(2,3,:) = M(1,3,:).*M(2,1,:)-M(1,1,:).*M(2,3,:);
    
    
    D(3,1,:) = M(2,1,:).*M(3,2,:)-M(2,2,:).*M(3,1,:);
    D(3,2,:) = M(1,2,:).*M(3,1,:)-M(1,1,:).*M(3,2,:);
    D(3,3,:) = M(1,1,:).*M(2,2,:)-M(1,2,:).*M(2,1,:);
    
    D = bsxfun(@rdivide,D,detM);
    D = reshape(D,[3 3 sM(3:end)]);
end

end