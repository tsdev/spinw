function oHandle = sw_getobject(tagName, fHandle)
% oHandle = SW_GETOBJECT(tagName, {fHandle}) returns all graphical objects
% those 'Tag' property string contains tagName string and is on figure
% fHandle.
%
% Input:
%
% fHandle       Handle of figure, optional. If undefined, objects on all
%               figures will be searched.
% tagName       String of the tag that the objects have, it can be a cell
%               that contains several strings.
%
% Output:
%
% oHandle       Handles of the found graphical objects in vector.
%
% See also SW, SW_GETFIGHANDLE, SW.PLOT.
%

if nargin == 0
    help sw_getobject;
    return
end

if iscell(tagName)
    oHandle = zeros(0,1);
    for ii = 1:numel(tagName)
        oHandle = [oHandle; findobj('-regexp','Tag',tagName{ii})]; %#ok<AGROW>
    end
else
    oHandle = findobj('-regexp','Tag',tagName);
end

if isempty(oHandle) || nargin == 1
    return;
end

hTemp = oHandle;
for ii = 1:length(oHandle)
    
    % if the object is a figure or figure descendent, return the figure.
    % Otherwise return 0.
    while any(hTemp(ii)) && ~strcmp('figure', get(hTemp(ii),'type'))
        hTemp(ii) = get(hTemp(ii),'parent');
    end
end

oHandle = oHandle(hTemp == fHandle);
end