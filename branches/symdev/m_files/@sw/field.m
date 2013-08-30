function varargout = field(obj,varargin)
% get/set magnetic field value
%
% {obj} = SW.FIELD(obj, B)
%
% If B is defined, it sets the magnetic field stored in sw object to B,
% where B is a 1x3 vector.
%
% B = SW.FIELD
%
% It returns the actual B field value.
%

if nargin == 1
    varargout{1} = obj.single_ion.field;
elseif nargin == 2
    B = varargin{1};
    if numel(B) == 3
        obj.single_ion.field = B(:)';
    else
        error('sw:magfield:ArraySize','Input magnetic field has to be a 3 element vector!');
    end
    if nargout > 0
        varargout{1} = obj;
    end
end

end % .magfield