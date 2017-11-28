function data = getdata(varargin)
% gets the data stored in an swplot figure
% 
% ### Syntax
% 
% `data = swplot.getdata`
% 
% `data = swplot.getdata(hFigure)`
%
% `data = swplot.getdata(field)`
%
% ### Description
% 
% `data = swplot.getdata` gets all the object data stored in the active
% swplot figure.
%
% `data = swplot.getdata(hFigure)` get all object data stored in the swplot
% figure identified by the `hFigure` handle.
%
% `data = swplot.getdata(field)` loads only the given field of the data
% structure.
% 
% ### Examples
%
% This example shows how the data of all objects on a 3D SpinW plot can be
% retrieved.
%
% ```
% >>model = sw_model('triAF',1)
% >>plot(model)
% >>swplot.getdata>>
% ```
% 
% ### Input Arguments
% 
% `hFigure`
% : Handle of the swplot figure, default value is the active figure.
% 
% `field`
% : String, determines the requested field name. If omitted, all
%   stored fields are returned.
% 
% ### See Also
% 
% [matlab.getappdata]
%

if nargin == 0
    hFigure = swplot.activefigure;
    arg = {hFigure};
elseif nargin == 1 && ischar(varargin{1})
    hFigure = swplot.activefigure;
    arg = {hFigure varargin{1}};
elseif nargin == 1
    arg = {gcf varargin{1}};
else
    arg = varargin;
end

data = getappdata(arg{:});

end