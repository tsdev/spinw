function hFigure = sw_getfighandle(tagLabel)
% hFigure = SW_GETFIGHANDLE(tagLabel) gives the handle of an active figure
% window that has the 'Tag' property identical to tagLabel. If possible, it
% returns the handle of the current figure. hFigure is empty if no figure
% with the necessary 'Tag' field is found.
%
% See also SW_ANNEALFIGURE.
%

if nargin == 0
    help sw_getfighandle;
    return
end

% List of figures that have tagLabel tag.
hList   = findobj('Tag',tagLabel);

if ~isempty(get(0,'Children')) && strcmp(tagLabel,get(gcf,'Tag'))
    hFigure = gcf;
elseif ~isempty(hList)
    hFigure = hList(end);
    figure(hFigure);
else
    hFigure = [];
end

end