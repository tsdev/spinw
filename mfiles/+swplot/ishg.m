function ishg = ishg(hFigure)
% does the swplot figure uses hgtransform
%
% ishg = SWPLOT.ISHG({hFigure})
%
% Input:
%
% hFigure       Handle of the swplot figure. Default is the selected
%               figure.
%

if nargin == 0
    hFigure = swplot.activefigure;
end

ishg = true(1,numel(hFigure));

for ii = 1:numel(hFigure)
    ishg(ii) = ~isempty(getappdata(hFigure(ii),'h'));
end

end