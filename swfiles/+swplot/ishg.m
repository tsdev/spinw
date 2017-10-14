function ishg = ishg(hFigure)
% does the swplot figure uses hgtransform
% 
% ### Syntax
% 
% `ishg = swplot.ishg`
%
% `ishg = swplot.ishg(hFigure)`
% 
% ### Description
% 
% `ishg = swplot.ishg` `true` if the active swplot figure uses
% [matlab.hgtransfrom], otherwise `false`.
% 
% `ishg = swplot.ishg(hFigure)` applies to the swplot figure referenced by
% the `hFigure` handle.
%
% ### See Also
%
% [swplot.figure]
%

if nargin == 0
    hFigure = swplot.activefigure;
end

ishg = true(1,numel(hFigure));

for ii = 1:numel(hFigure)
    ishg(ii) = ~isempty(getappdata(hFigure(ii),'h'));
end

end