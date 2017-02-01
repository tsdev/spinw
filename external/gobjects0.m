function out = gobjects(varargin)
% use double insted of graphics handle class for Matlab prior 8.1

if verLessThan('matlab','8.1')
    out = zeros(varargin{:});
else
    out = gobjects(varargin{:});
end

end