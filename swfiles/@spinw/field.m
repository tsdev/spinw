function varargout = field(obj,varargin)
% get/set magnetic field value
% 
% ### Syntax
% 
% `field(obj,B)`
% `B = field(obj)`
% 
% ### Description
% 
% `field(obj,B)` sets the magnetic field stored in `obj.single_ion.field`
% to `B`, where `B` is a $[1\times 3]$ vector.
%  
% `B = field(obj)` returns the current value of the magnetic field value
% stored in `obj`.
%  
% ### See Also
% 
% [spinw] \| [spinw.temperature] \| [spinw.single_ion]
%

if nargin == 1
    varargout{1} = obj.single_ion.field;
elseif nargin == 2
    B = varargin{1};
    if numel(B) == 3
        if obj.symbolic
            if isa(B,'sym')
                obj.single_ion.field = B(:)';
            else
                obj.single_ion.field = B(:)'*sym('B','real');
            end
        else
            obj.single_ion.field = B(:)';
        end
    else
        error('spinw:magfield:ArraySize','Input magnetic field has to be a 3 element vector!');
    end
    if nargout > 0
        varargout{1} = obj;
    end
end

end % .field