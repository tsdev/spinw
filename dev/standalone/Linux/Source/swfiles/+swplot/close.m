function close(varargin)
% close swplot figure
%
% SWPLOT.CLOSE()
%
% Closes the active swplot figure.
%
% SWPLOT.CLOSE('all')
%
% Closes all swplot figure.
%
% SWPLOT.CLOSE(hFigure)
%
% closes the swplot figure of the given handle (hFigure).
%
% Input:
%
% hFigure   Handle of the swplot figure window.
%
% See also SWPLOT.FIGURE.
%

if nargin == 0
    % check if there is any swplot figure
    activeTag = swpref.getpref('tag',[]);
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
    activeTag = swpref.getpref('tag',[]);
    % tag for inactive figures
    inactiveTag = ['inactive_' activeTag];
    
    % find and close all swplot figures
    close(findobj('tag',activeTag));
    close(findobj('tag',inactiveTag));
end

end