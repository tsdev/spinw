function varargout = notwin(obj)
% removes any twin added to the sw object
%
% NOTWIN(obj)
%
% The function keeps only the original twin.
%

obj.twin.vol = 1;
obj.twin.rotc = eye(3);

if nargout >0
    varargout{1} = obj;
end

end