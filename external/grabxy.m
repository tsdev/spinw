function pts = grabxy(fName, ax, logax)
% reads coordinates from raster image
%
% pts = GRABXY(fName, ax, {logax})
%
% Input:
%
% fName         String, path to the image file.
% ax            (x,y) coordinates of the three axis points, with dimensions
%               of 2x3. If ax is omitted or empty, GRABXY just shows the
%               image.
% logax         Optional input string, cn be used if the axis is in
%               logaritmic units, possible options:
%                   'logx'
%                   'logy'
%                   'logxy'
%
% Output:
%
% pts           Matrix, where the first row are the x coordinates the
%               second row is the y coordinates.
%
% How to use the function?
%
% You need an image where you know the coordinates of three points. These
% coordinates has to be given as the ax input to the grabxy() function.
% After calling the function, it shows the image file and the first three
% mouse clicks on the figure will determine the position of the three
% predefined points (shown as red circles) given in ax. Any other click
% afterwards will grab coordinates (shown as green circles). The coordinates
% of the green circles will be returned. If you don't want to grab more
% points, just right click anywhere on the figure.
%

if nargin == 0
    help grabxy
    return
end

if nargin < 3
    logax = '';
end

% open image file
img = imread(fName);

% display the image
hFig = figure('ToolBar','none','MenuBar','none','name','GrabXY');
image(img);
axis off

% just show the image if no ax is provided
if nargin == 1 || (nargin==2 && isempty(ax))
    while 1
        [~, ~, but] = ginput(1);      % get a point
        if ~isequal(but, 1)             % stop if not button 1
            break
        end
    end
    close(hFig);
    if nargout > 0
        pts = zeros(2,0);
    end
    return
end

idx = 0;
hold on
pts0 = zeros(2,0);

while 1
    [xi, yi, but] = ginput(1);      % get a point
    if ~isequal(but, 1)             % stop if not button 1
        break
    end
    idx = idx + 1;
    pts0(1,idx) = xi;
    pts0(2,idx) = yi;
    
    if idx < 4
        col = 'r';
    else
        col = 'g';
    end
    
    if idx>4
        % draw line
        plot(pts0(1,[-1 0]+idx), pts0(2,[-1 0]+idx), [col 'o-']);
    else
        % only pints
        plot(pts0(1,idx), pts0(2,idx), [col 'o']);
    end
end

hold off

% solve linear equation
if any(size(ax)-[2 3])
    if any(size(ax)-[3 2])
        error('Wrong dimensions of ax, has to be either 2x3 or 3x2!');
    else
        ax = ax';
    end
end

% linear equation: B = A*X
% define matrices
B = [ax(1,:) ax(2,:)]';

if idx>4
    As = [pts0(:,1)' 1;pts0(:,2)' 1;pts0(:,3)' 1];
    A = zeros(6);
    A(1:3,1:3) = As;
    A(4:6,4:6) = As;
    
    % solve linear equations
    X = A\B;
    
    % tranformation matrix x'=Tx + x0
    T = [X(1:2)';X(4:5)'];
    x0 = X([3 6]);
    
    % new transformed coordinates
    pts = bsxfun(@plus,T*pts0(:,4:end),x0);
else
    pts = zeros(2,0);
end

% convert values into logaritmic units
switch logax
    case 'logx'
        idx = 1;
    case 'logy'
        idx = 2;
    case 'logxy'
        idx = [1 2];
    case ''
        idx = [];
end

for ii = 1:numel(idx)
        idx0 = idx(ii);
        coo = pts(idx0,:);
        % limiting points
        aLim = [min(ax(idx0,:)) max(ax(idx0,:))];
        A = aLim(1)-aLim(2)/log(aLim(1)/aLim(2));
        B = aLim(1)-A*log(aLim(1));
        pts(idx0,:) = exp((coo-B)/A);
end

% close figure
close(hFig)

end