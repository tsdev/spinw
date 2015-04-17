function varargout = findmin(X,varargin)
% find the minimum valued element in a vector
%
% I = FINDMIN(X,val)
% I = FINDMIN(X,val,k)
% I = FINDMIN(X,val,k,direction)
%
% [row,col] = FINDMIN(___)
% [row,col,v] = FINDMIN(___)
%
% It has the same parameters and output as find().
%
% See also find.

switch nargout
    case 0
        varargout{1} = find(X==min(X),varargin{:});
    case 1
        varargout{1} = find(X==min(X),varargin{:});
    case 2
        [varargout{1},varargout{2}] = find(X==min(X),varargin{:});
    case 3
        [row,col] = find(X==min(X),varargin{:});
        
        varargout{1} = row;
        varargout{2} = col;
        varargout{3} = X(sub2ind(size(X),row,col));
end

end