function clear(varargin)
% clear swplot figure
%
% SWPLOT.CLEAR(hFigure)
%
% clears swplot figure correspondign to hFigure handle
%
% SWPLOT.CLEAR
%
% clears the active swplot figure

swplot.delete(varargin{:},0)

end