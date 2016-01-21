function [R2, km] = ksearch(abc, k, N, kR, kD)
% search for best magnetic k-vector
%
% [km, R2] = KSEARCH(abc, km, Ngrid, kmRange, kmNum)
%
% Input:
%
% abc   Vector defining the crystal lattice with the following numbers:
%           [a b c alpha beta gamma], angles are in degree. If the crystal
%           is defined in SpinW, just use the obj.abc command.
% km    List of magnetic Bragg peak positions in A-1.
% Ngrid Vector with 3 numbers or scalar, the number of reciprocal lattice
%       points along H, K and L to test. If a single number is given, the
%       same number of points will be tested along all three dimensions.
%       The HKL grid will be: -Ngrid:Ngrid.
% kmRange   Either a vector with 3 elements or a scalar, defines the range
%       of magnetic ordering wave vectors to be tested. For example:
%           kmRange = 0.3, then ordering wave vectors will be search
%           between 0 and 0.3.
% kmNum The number of ordering wave vectors between 0 and kmRange.
%
% Output:
%
% km    Matrix with dimensions 10x3, returns the best 10 magnetic ordering
%       wave vectors.
% R2    Column vector containing the R^2 values for the best 10 fits.
%


if numel(N) == 1
    N = [N N N];
end

if numel(kR) == 1
    kR = [kR kR kR];
end

if numel(kD) == 1
    kD = [kD kD kD];
end

km = cell(1,3);

for ii = 1:3
    km{ii} = [linspace(0,kR(ii),kD(ii)) 1/2 1/3];
end

[km{1},km{2},km{3}] = ndgrid(km{:});
km = reshape(cat(4,km{:}),[],3);

obj = sw;
obj.genlattice('lat_const',abc(1:3),'angled',abc(4:6));

T  = 2*pi*inv(obj.basisvector); %#ok<MINV>
R2 = zeros(size(km,1),1);

k  = k(:);
nk = numel(k);

H = (-N(1):N(1));
K = (-N(2):N(2));
L = (-N(3):N(3));

hkl0 = cell(1,3);
[hkl0{1}, hkl0{2},hkl0{3}] = ndgrid(H,K,L);

hkl0 = cat(4,hkl0{:});

for jj = 1:size(km,1)
    
    hkl = [bsxfun(@minus,hkl0,permute(km(jj,:),[1 3 4 2])); bsxfun(@plus,hkl0,permute(km(jj,:),[1 3 4 2]))];
    
    hklA = sqrt(sum(sum(bsxfun(@times,hkl,permute(T,[3 4 5 1 2])),4).^2,5));
    
    dQ = reshape((bsxfun(@minus,hklA,permute(k,[2 3 4 1]))).^2,[],nk);
    
    R2(jj) = sum(min(dQ,[],1));
    
end

% find the best solutions
[R2,idx] = sort(R2);
km = km(idx,:);

R2 = R2(1:10)*1e4;
km = km(1:10,:);

end