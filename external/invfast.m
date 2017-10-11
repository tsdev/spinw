function I = invfast(M,opt)
% fast calculation of the matrix inverse of stacked matrices
%
% I = INVFAST(M,{opt})
%
% The function is vectorised and efficient for very large number of small
% matrices (1x1, 2x2 and 3x3). For larger matrices the code calls the
% built-in inv() function. If opt is 'diag', only the diagonal of the
% inverse matrices are calculated.
%
% Input:
%
% M         Matrix with dimensions of [D,D,N1,N2,...], where D can be
%           1, 2 or 3.
% opt       String, with possile values:
%               'full'  Calculate the inverse matrices, the output matrix
%                       will have the same dimensions as the input matrix.
%                       Default option.
%               'diag'  Calculate only the diagonal elements of the inverse
%                       matrix, the output matrix will have dimensions of
%                       [D,N1,N2,...].
%               'sum'   Calculate only the sum of the elements of the
%                       inverse matrix, the output will have dimensions of
%                       [1,N1,N2,...].
%           Default value is 'full'.
%
% Output:
%
% I         Matrix with dimensions determined by opt.
%

if nargin == 0
    help invfast
    return
end

if nargin < 2
    opt = 'full';
end

sM = [size(M) 1];

% matrix dimension
D = sM(1);

if sM(1)~= sM(2)
    error('invfast:WrongInput','Input matrix has wrong dimensions!')
end

% deal with vectors
M    = reshape(M,D,D,[]);
nMat = size(M,3);

if D>3
    % call inv for large matrices, no speedup
    switch opt
        case 'diag'
            I = zeros(D,nMat);
            for ii = 1:nMat
                T       = inv(M(:,:,ii));
                I(:,ii) = diag(T);
            end
            I = reshape(I,[D sM(3:end)]);
        case 'full'
            I = zeros(size(M));
            for ii = 1:nMat
                I(:,:,ii) = inv(M(:,:,ii));
            end
            I = reshape(I,[D D sM(3:end)]);
        case 'sum'
            I = zeros(1,nMat);
            for ii = 1:nMat
                T       = inv(M(:,:,ii));
                I(1,ii) = sum(T(:));
            end
            I = reshape(I,[1 sM(3:end)]);
    end
    return
end

% calculation of the matrix determinant
switch D
    case 2
        detM = M(1,1,:).*M(2,2,:)-M(1,2,:).*M(2,1,:);
    case 3
        detM = M(1,1,:).*(M(2,2,:).*M(3,3,:)-M(2,3,:).*M(3,2,:))...
            +M(1,2,:).*(M(2,3,:).*M(3,1,:)-M(2,1,:).*M(3,3,:))...
            +M(1,3,:).*(M(2,1,:).*M(3,2,:)-M(2,2,:).*M(3,1,:));
end

switch opt
    case 'diag'
        % calculate only the diagonal elements
        switch D
            case 1
                I = 1./M;
                
            case 2
                I = zeros(2,nMat);
                I(1,:) = M(2,2,:);
                I(2,:) = M(1,1,:);
                I = bsxfun(@rdivide,I,permute(detM,[1 3 2]));
                
            case 3
                I = zeros(3,nMat);
                I(1,:) = M(2,2,:).*M(3,3,:)-M(2,3,:).*M(3,2,:);
                I(2,:) = M(1,1,:).*M(3,3,:)-M(1,3,:).*M(3,1,:);
                I(3,:) = M(1,1,:).*M(2,2,:)-M(1,2,:).*M(2,1,:);
                I = bsxfun(@rdivide,I,permute(detM,[1 3 2]));
                
        end
        I = reshape(I,[D sM(3:end)]);
    case 'full'
        % calculate vectorized inverse
        switch D
            case 1
                I = 1./M;
            case 2
                I = zeros(2,2,nMat);
                I(1,1,:) =  M(2,2,:);
                I(1,2,:) = -M(1,2,:);
                I(2,1,:) = -M(2,1,:);
                I(2,2,:) =  M(1,1,:);
                I = bsxfun(@rdivide,I,detM);
                
            case 3
                I = zeros(3,3,nMat);
                
                I(1,1,:) = M(2,2,:).*M(3,3,:)-M(2,3,:).*M(3,2,:);
                I(1,2,:) = M(1,3,:).*M(3,2,:)-M(1,2,:).*M(3,3,:);
                I(1,3,:) = M(1,2,:).*M(2,3,:)-M(1,3,:).*M(2,2,:);
                
                
                I(2,1,:) = M(2,3,:).*M(3,1,:)-M(2,1,:).*M(3,3,:);
                I(2,2,:) = M(1,1,:).*M(3,3,:)-M(1,3,:).*M(3,1,:);
                I(2,3,:) = M(1,3,:).*M(2,1,:)-M(1,1,:).*M(2,3,:);
                
                
                I(3,1,:) = M(2,1,:).*M(3,2,:)-M(2,2,:).*M(3,1,:);
                I(3,2,:) = M(1,2,:).*M(3,1,:)-M(1,1,:).*M(3,2,:);
                I(3,3,:) = M(1,1,:).*M(2,2,:)-M(1,2,:).*M(2,1,:);
                
                I = bsxfun(@rdivide,I,detM);
        end
        I = reshape(I,[D D sM(3:end)]);
    case 'sum'
        % calculate the sum of the inverse
        switch D
            case 1
                I = 1./M;
            case 2
                I =  M(1,1,:)+M(2,2,:)-M(1,2,:)-M(2,1,:);
                I = permute(bsxfun(@rdivide,I,detM),[1 3 2]);
                
            case 3
                
                I = M(2,2,:).*M(3,3,:)-M(2,3,:).*M(3,2,:)+ M(1,3,:).*M(3,2,:)-M(1,2,:).*M(3,3,:) + ...
                    M(1,2,:).*M(2,3,:)-M(1,3,:).*M(2,2,:) + M(2,3,:).*M(3,1,:)-M(2,1,:).*M(3,3,:) + ...
                    M(1,1,:).*M(3,3,:)-M(1,3,:).*M(3,1,:) + M(1,3,:).*M(2,1,:)-M(1,1,:).*M(2,3,:) + ...
                    M(2,1,:).*M(3,2,:)-M(2,2,:).*M(3,1,:) + M(1,2,:).*M(3,1,:)-M(1,1,:).*M(3,2,:) + ...
                    M(1,1,:).*M(2,2,:)-M(1,2,:).*M(2,1,:);
                
                I = permute(bsxfun(@rdivide,I,detM),[1 3 2]);
        end
        I = reshape(I,[1 sM(3:end)]);
end

end