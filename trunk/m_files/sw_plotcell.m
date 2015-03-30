function sw_plotcell(dat,xLab,yLab,tLab,aRange)
% plots cell structure with circles


nDat = numel(dat.x);
zMax = -inf;

% find maximum data value
for ii = 1:nDat
    zMax = max([abs(real(dat.z{ii}(:))); zMax]);
end

if nargin>4
    xMin = aRange(1);
    xMax = aRange(2);
    yMin = aRange(3);
    yMax = aRange(4);
else
    xMin = min(dat.x)-max(diff(dat.x));
    xMax = max(dat.x)+max(diff(dat.x));
    yMin = inf;
    yMax = -inf;
    
    % find x and y range
    for ii = 1:nDat
        yMin = min([abs(real(dat.y{ii}(:))); yMin]);
        yMax = max([abs(real(dat.y{ii}(:))); yMax]);
    end
end

% radius of circles
rx = (xMax-xMin)/50;
ry = (yMax-yMin)/50;

C = sw_circle([0 0 0]',[0 0 1]',1,300);
Cx = C(1,:)*rx/zMax;
Cy = C(2,:)*ry/zMax;

% draw circles
for ii = 1:nDat
    for jj = 1:numel(dat.y{ii})
        r = dat.z{ii}(jj);
        p = patch(dat.x(ii)+Cx*r,dat.y{ii}(jj)+Cy*r,'k');
        hold on
    end
end

dy = (yMax-yMin)/20;
axis([xMin xMax yMin-dy yMax+dy]);

if nargin>1
    xlabel(xLab)
end
if nargin>2
    ylabel(yLab)
end
if nargin>3
    title(tLab)
end

end