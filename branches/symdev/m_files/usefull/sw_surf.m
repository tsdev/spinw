function sw_surf(X,Y,C,cLim,cMap,maxPatch)
% SW_SURF draws a two dimensional plot using the patch function when the
% figure is saved in .pdf format, it is saved as a vector image
%
% Usage:
%
% SW_SURF(C)
%
% SW_SURF(X,Y,C,cLim,cMap)
%
% Input:
%
% X         Matrix of x coordinates, dimensions are [nX nY].
% Y         Matrix of y coordinates, dimensions are [nX nY].
% C         Matrix to plot, dimensions are [nX nY].
% cLim      Color axis limits, if undefined the minimum and maximum values
%           C are taken as limits.
% cMap      Colormap, default is @jet.
% maxPatch  Maximum number of pixels that can be plotted using the patch()
%           function. Using patch for color plot can be slow on older
%           machines, but the figure can be exported afterwards as a vector
%           graphics, using the print() function. Default is 1000.
%
% Change of the color axis after plotting is not possible, replot is
% necessary.
%

if nargin < 6
    maxPatch = 1000;
end

if nargin == 1
    C = X;
    xAxis = 1:size(C,1);
    yAxis = 1:size(C,2);
    [X, Y] = ndgrid(xAxis,yAxis);
end

if nargin > 4
    N = numel(cMap)/3;
else
    N = 64;
    cMap = jet(N);
end

if nargin < 4
    
    idxD = C-min(C(:));
    
    idxD = idxD/max(abs(idxD(:)))*(N-1)+1;
    cLim = [min(C(:)) max(C(:))];
else
    idxD = (C-cLim(1))/(cLim(2)-cLim(1))*(N-1)+1;
end

if numel(C) < maxPatch
    idxD = floor(idxD);
    idxD(idxD<1) = 1;
    idxD(idxD>N) = N;
    
    dx = diff(X,1,2);
    dx = [dx(:,1) dx];
    dy = diff(Y);
    dy = [dy(1,:);  dy];
    
    for ii = 1:size(C,1)
        for jj = 1:size(C,2)
            patch([-1 1 1 -1]*dx(ii,jj)/2+X(ii,jj),[-1 -1 1 1]*dy(ii,jj)/2+Y(ii,jj),...
                cMap(floor(idxD(ii,jj)),:),'edgecolor',cMap(idxD(ii,jj),:));
        end
    end
    
else
    hSurf = surf(X,Y,C);
    view(2)
    set(hSurf,'EdgeAlpha',0);
    colormap(cMap);
    caxis(cLim)
end

end