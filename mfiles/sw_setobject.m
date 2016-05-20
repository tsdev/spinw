function varargout = sw_getobject(tagName, varargin)
% sets graphic object properties that has the given 'Tag' property
%
% oHandle = SW_SETOBJECT(tagName, pName, pValue, ..., {fHandle}) 
%
% The function sets properties of graphical objects those 'Tag' property
% string contains tagName string and is on figure fHandle.
%
% Input:
%
% tagName       String of the tag that the objects have, it can be a cell
%               that contains several strings.
% fHandle       Handle of figure, optional. If undefined, objects on all
%               figures will be searched.
%
% Output:
%
% oHandle       Handles of the found graphical objects in vector. Optional.
%
% See also SW, SW_GETFIGHANDLE, SW.PLOT, SW_GETOBJECT.
%

if nargin == 0
    help sw_setobject
    return
end

if mod(nargin,2) == 0
    % last argument is figure handle
    oHandle = sw_getobject(tagName,varargin{end});
    varargin = varargin(1:(end-1));
else
    oHandle = sw_getobject(tagName);
end

set(oHandle,varargin{:});

if nargout == 1
    varargout{1} = oHandle;
end

end