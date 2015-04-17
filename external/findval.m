function varargout = findval(X,val,varargin)
% find the an element in a vector that is closest to a given value
%
% I = FINDVAL(X,val)
% I = FINDVAL(X,val,k)
% I = FINDVAL(X,val,k,direction)
%
% [row,col] = FINDVAL(___)
% [row,col,v] = FINDVAL(___)
%
% It has the same parameters and output as find() plus the extra val value.
%
% See also find.

switch nargout
    case 0
        varargout{1} = find(abs(X-val)==min(abs(X-val)),varargin{:});
    case 1
        varargout{1} = find(abs(X-val)==min(abs(X-val)),varargin{:});
    case 2
        [varargout{1},varargout{2}] = find(abs(X-val)==min(abs(X-val)),varargin{:});
    case 3
        [row,col] = find(abs(X-val)==min(abs(X-val)),varargin{:});
        
        varargout{1} = row;
        varargout{2} = col;
        varargout{3} = X(sub2ind(size(X),row,col));
end

end