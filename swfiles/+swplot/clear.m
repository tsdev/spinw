function clear(varargin)
% clears swplot figure
% 
% ### Syntax
% 
% `swplot.clear`
%
% `swplot.clear(hFigure)`
% 
% ### Description
% 
% `swplot.clear` clears the active swplot figure.
%
% `swplot.clear(hFigure)` clears the swplot figure correspondign to
% `hFigure` handle
%  
% ### See Also
% 
% [swplot.figure]
%

swplot.delete(varargin{:},0)

end