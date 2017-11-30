function close(varargin)
% closes swplot figure
% 
% ### Syntax
% 
% `swplot.close`
% 
% `swplot.close('all')`
%
% `swplot.close(hFigure)`
% 
% ### Description
% 
% `swplot.close` closes the active swplot figure.
%  
% `swplot.close('all')` closes all swplot figure.
%  
% `swplot.close(hFigure)` closes the swplot figure corresponding to
% `hFigure` handle.
% 
% ### See Also
% 
% [swplot.figure]
%

pref = swpref;


if nargin == 0
    % check if there is any swplot figure
    activeTag = pref.tag;
    inactiveTag = ['inactive_' activeTag];
    if isempty([findobj('tag',activeTag) findobj('tag',inactiveTag)])
        % nothing to close
        return
    end
    
    hFigure = swplot.activefigure;
else
    if strcmp('all',varargin{1})
        hFigure = [];
    else
        hFigure = varargin{1};
    end
end

if ~isempty(hFigure)
    close(hFigure);
else
    % close all swplot figure
    activeTag = pref.tag;
    % tag for inactive figures
    inactiveTag = ['inactive_' activeTag];
    
    % find and close all swplot figures
    close(findobj('tag',activeTag));
    close(findobj('tag',inactiveTag));
end

end