function A = sumn(A,dim, varargin)
% sum of array elements along multiple dimensions
%
% S = SUMN(A,dim,...)
%
% Works the same way as the Matlab built-in sum() function, but it can do
% summation along multiple dimensions with a single call.
%
% See also SUM.
%

if nargin == 0
    help sumn
    return
end

if nargin == 1
    % sum along the first non-singleton dimension
    dim = find(size(A)>1,1,'first');
end

% loop over all summed dimensions
for ii = 1:numel(dim)
    A = sum(A,dim(ii),varargin{:});
end

end