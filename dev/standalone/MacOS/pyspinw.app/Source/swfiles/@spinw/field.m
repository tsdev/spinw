function varargout = field(obj,varargin)
% get/set magnetic field value
%
% FIELD(obj, B)
%
% If B is defined, it sets the magnetic field stored in spinw object to B,
% where B is a 1x3 vector.
%
% B = FIELD(obj)
%
% The function returns the current B field value stored in obj.
%
% See also SPINW, SPINW.TEMPERATURE.
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