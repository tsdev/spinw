function varargout = temperature(obj,varargin)
% get/set stored temperature value
%
% TEMPERATURE(obj, T)
%
% If T is defined, it sets the temperature stored in obj object
% to T, where T is scalar. The units of temerature is
% determined by the spinw.unit.kB value, default is Kelvin.
%
% T = TEMPERATURE(obj)
%
% The function returns the current temperature value stored in
% obj.
%

if nargin == 1
    varargout{1} = obj.single_ion.T;
elseif nargin == 2
    T = varargin{1};
    if numel(T) == 1
        obj.single_ion.T = T;
    else
        error('spinw:temperature:ArraySize','Input temperature has to be scalar!');
    end
    if nargout > 0
        varargout{1} = obj;
    end
end

end