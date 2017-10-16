function varargout = subplot(varargin)
% create subplots with variable gaps between axes
% 
% ### Syntax
% 
% `swplot.subplot(m,n,p,space)`
% 
% `swplot.subplot([m,n,p],space)`
%
% ### Description
% 
% `swplot.subplot(m,n,p,space)` is equivalent to the [matlab.subplot]
% command, except that the space between axes can be controlled.
% 
% ### Input Arguments
% 
% `m,n,p`
% : Three integer that defines subplot, for details see the
%   built-in [matlab.subplot] command.
% 
% `space`
% : Vector with elements: `[margin hgap vgap]`, where:
%   * `margin`  Top and right margin at the figure edge.
%   * `hgap`    Left margin and horizontal gap between axes.
%   * `vgap`    Bottom margin and vertical gap between axes.
% 
% ### See Also
% 
% [matlab.subplot]
%

if nargin==3 || nargin == 1
    space = [0 0 0];
else
    space = varargin{end};
end

if nargin<3
    m = varargin{1}(1);
    n = varargin{1}(2);
    p = varargin{1}(3);
else
    m = varargin{1};
    n = varargin{2};
    p = varargin{3};
end

if p>m*n || p<0
    error('subplot:WrongInput','Subpanel index is out of range!')
end

% space = [margin hgap vgap]
margin = space(1);
hgap   = space(2);
vgap   = space(3);

posi = [mod(p-1,n) floor((p-1)/n)];

w = ((1-margin-hgap)-(n-1)*hgap)/n;
h = ((1-margin-vgap)-(m-1)*vgap)/m;

pos = [hgap+posi(1)*(hgap+w) 1-(margin+(posi(2)+1)*h+posi(2)*vgap) w h];

varargout{:} = subplot('Position',pos);

end