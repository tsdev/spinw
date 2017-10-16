function subfigure(m,n,p,hFigure)
% changes position of figure window on the screen
% 
% ### Syntax
% 
% `swplot.subfigure(m,n,p)`
% 
% `swplot.subfigure(m,n,p,hFigure)`
%
% ### Description
% 
% `swplot.subfigure(m,n,p)` changes the position of the current figure
% window on the screen, the position is determined similarly to the Matlab
% function [matlab.subplot]. Here the screen is the canvas where the figure
% window is positioned.
%  
% The function divides the display into an $m$-by-$n$ grid and moves the
% figure window in the position specified by $p$. It numbers the figures by
% row major, such that the first figure is the first column of the first
% row, the second figure is the second column of the first row, and so on.
%
% `swplot.subfigure(m,n,p,hFigure)` repositions the figure related to
% `hFigure` handle.
% 
% ### Input Arguments
% 
% `m,n,p`
% : Integer numbers that define the figure window position.
% 
% `hFigure`
% : Handle of the figure window, optional. Default value is [matlab.gcf].
%

if nargin < 4
    % find active figure
    hFigure = gcf;
end

position = [m n p];

if p<1 || p>m*n
    error('subfigure:WrongInput','Value of p is illegal!');
end

% Position figure window on the screen
fUnit = get(hFigure,'Units');
set(hFigure,'Units','pixels');

idx = position(3);
nx  = position(2);
ny  = position(1);
position = [ceil(idx/nx) mod(idx-1,nx)+1 ny nx];

if numel(position)==4
    % get screen size
    unit0 = get(0,'units');
    set(0,'units','pixels');
    scSize = get(0,'screensize');
    set(0,'units',unit0);
    
    fWidth  = scSize(3)/position(4);
    fHeight = scSize(4)/position(3);
else
    fPos = get(hFigure,'outerPosition');
    fWidth  = fPos(3);
    fHeight = fPos(4);
end

% Display size in pixel
dSize = get(0,'ScreenSize');
set(hFigure,'outerPosition',[(position(2)-1)*fWidth dSize(4)-20-position(1)*fHeight fWidth fHeight]);
set(hFigure,'Units',fUnit);

end