function S = sumargs(varargin)
% sums up all arguments

S = varargin{1};

for ii = 2:nargin
    S = S + varargin{ii};
end

end