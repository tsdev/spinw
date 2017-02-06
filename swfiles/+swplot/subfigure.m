function subfigure(m,n,p,hFigure)
% changes position of figure window on the screen
%
% SWPLOT.SUBFIGURE(m,n,p,{hFigure})
%
% Input:
%
% Changes the position of the figure window on the screen, the position is
% determined similarly to the Matlab function subplot(). Here the screen is
% the canvas where the figure window is positioned.
%
% SUBFIGURE divides the display into an m-by-n grid and moves the figure
% window in the position specified by p. It numbers its figures by row,
% such that the first figure is the first column of the first row, the
% second figure is the second column of the first row, and so on.
%
% Input:
%
% m,n,p     Integer numbers, defining figure window position.
% hFigure   Handle of the figure window, optional. Default is gcf.
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