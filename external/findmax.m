function varargout = findmax(X,varargin)
% find the maximum valued element in a vector
%
% I = FINDMAX(X,val)
% I = FINDMAX(X,val,k)
% I = FINDMAX(X,val,k,direction)
%
% [row,col] = FINDMAX(___)
% [row,col,v] = FINDMAX(___)
%
% It has the same parameters and output as find().
%
% See also find.

switch nargout
    case 0
        varargout{1} = find(X==max(X),varargin{:});
    case 1
        varargout{1} = find(X==max(X),varargin{:});
    case 2
        [varargout{1},varargout{2}] = find(X==max(X),varargin{:});
    case 3
        [row,col] = find(X==max(X),varargin{:});
        
        varargout{1} = row;
        varargout{2} = col;
        varargout{3} = X(sub2ind(size(X),row,col));
end

end
