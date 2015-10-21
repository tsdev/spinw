function subwin(h,w,p)
% positions of the figure window on the screen.
%
% subwin(m,n,p)
%
% Works similarly as subplot(), just positions figure windows relative to
% the screen instead of axes relative to the figure. Three number vector,
% [m n p]. It divides the screen into an m-by-n grid and positions the
% figure in the position specified by p. The window positions are numbered
% by row, such that the first window position is the first column of the
% first row, the second window position is the second column of the first
% row, and so on. If m or n equal to zero, the original size of the figure
% window is used to determine the grid.
%
% See also subplot().
%

fUnit = get(gcf,'Units');
set(gcf,'Units','pixels');

unit0 = get(0,'units');
set(0,'units','pixels');
scSize = get(0,'screensize');
set(0,'units',unit0);

% display size
dHeight = scSize(4);
dWidth  = scSize(3);

% determine the upper menu size
tFig = figure;
set(tFig,'units','pixels')
set(tFig,'outerposition',[1 dHeight 50 50])
drawnow
oP = get(tFig,'outerPosition');
close(tFig);

dHeight = oP(2)+oP(4);
if h==0 || w==0
    % use the size of figure
    fPos = get(gcf,'outerPosition');
    fHeight = fPos(4);
    fWidth  = fPos(3);
else
    % determine the size of figure from the display size
    fHeight = dHeight/h;
    fWidth  = dWidth/w;
end

% determine the number of window position on the screen
h0 = floor(dHeight/fHeight);
w0 = floor(dWidth/fWidth);

% position
px = mod(p-1,w0);
py = floor((p-1)/w0);

set(gcf,'outerPosition',[px*fWidth+1 (h0-py-1)*fHeight+dHeight-fHeight*h0 fWidth fHeight]);


set(gcf,'Units',fUnit);

end