function data = getdata(varargin)
% gets the data stored in an swplot figure
%
% data = SWPLOT.GETDATA({hFigure},{field}
%
% Input:
%
% hFigure       Handle of the swplot figure. Default is the active figure.
% field         String, determined the stored field name. If omitted, all
%               stored data are returned.
%
% See also GETAPPDATA.
%

if nargin == 0
    hFigure = swplot.activefigure;
    arg = {hFigure};
elseif nargin == 1 && isstring(varargin{1})
    hFigure = swplot.activefigure;
    arg = {hFigure varargin{1}};
elseif nargin == 1
    arg = {gcf varargin{1}};
else
    arg = varargin;
end

data = getappdata(arg{:});

end