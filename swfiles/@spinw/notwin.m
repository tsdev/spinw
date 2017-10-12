function varargout = notwin(obj)
% removes all twins
% 
% ### Syntax
% 
% `notwin(obj)`
% 
% ### Description
% 
% `notwin(obj)` removes any crystallographic twin added using the
% [spinw.addtwin] function.
%  
% ### See Also
%
% [spinw.addtwin]
%

obj.twin.vol = 1;
obj.twin.rotc = eye(3);

if nargout >0
    varargout{1} = obj;
end

end