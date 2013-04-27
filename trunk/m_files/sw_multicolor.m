function cMat = sw_multicolor(vMat, cMap, cLim, nCol, pflipud)
% cMat = sw_multicolor(vMat, cMap, cLim, {nCol}, {flipud}) creates RGB color data for
% multiple 2D overlapping plots.
%
% Input:
%
% vMat      Matrix that contains the input 2D data, dimensions are
%           [d1 d2 nPlot]. Where each plot has a dimensions of [d1 d2].
% cMap      Cell of colormap functions containing used for the different
%           overlayed plots. For example:
%           cMap = {@copper @gray}.
% cLim      Maximum and minimum values of the color map. Values in vMat
%           smaller than the minimum and larger than the maximum will be
%           shown with the minimum and maximum values in the colormap
%           respectively. The dimensions of cLim is [1 2].
% nCol      Number of colors in the colormap. Optional, default value is
%           100.
% flipud    If true the colormaps are inverted. Optional, default is false.
%
% Output:
%
% cMat      Matrix with equal dimensions to the input times three for the
%           red, green and blue channels, dimensions are [d1 d2 3].
%
% Example:
%           Plotting of two random matrices (dimensions are [100 100]) with
%           red and blue colors:
%               cMat = sw_multicolor(rand(100,100,2),[1 0;0 1;0 0],[0 1]);
%               image(cMat);
%

if nargin < 3
    help sw_multicolor;
    return;
end

if nargin == 3
    nCol = 100;
end
if nargin < 5
    pflipud = false;
end

nPlot = size(vMat,3);
d12   = size(vMat);
d12   = d12(1:2);

% scale values into the (0,1) range
vMat(vMat<cLim(1)) = cLim(1);
vMat(vMat>cLim(2)) = cLim(2);
vMat = (vMat-cLim(1))/(cLim(2)-cLim(1));
vMat = round(vMat*(nCol-1))+1;

% create colormaps for every individual plot
cMapGen = zeros(nCol,4,nPlot);
for ii = 1:nPlot
    if pflipud
        cMapGen(:,1:3,ii) = flipud(cMap{ii}(nCol));
    else
        cMapGen(:,1:3,ii) = cMap{ii}(nCol);
    end
end

% convert colormaps to CMYK colors
for ii = 1:nPlot
    cMapGen(:,:,ii) = rgb2cmyk(cMapGen(:,1:3,ii));
end

% get CMYK colors for every individual plot
vCMYK = zeros(d12(1),d12(2),4,nPlot);
for ii = 1:nPlot
    vCMYK(:,:,:,ii) = reshape(cMapGen(vMat(:,:,ii),:,ii),[d12 4]);
end

% mix CMYK colors
cMat = zeros(d12(1),d12(2),4);
for ii = 1:nPlot
    cMat = max(cMat,vCMYK(:,:,:,ii));
end

% convert CMYK colors back to RGB
cMat = reshape(rgb2cmyk(reshape(cMat,[],4)),[d12 3]);

end

function c = rgb2cmyk(c)
% c = RGB2CMYK(c)   Converts between RGB- and CMYK- Colors
%
% RedGreenBlue <---> CyanMangentaYellowBlack
%
% CMYK = RGB2CMYK(  RGB )
%
%  RGB = RGB2CMYK( CMYK )
%
% RGB :  N by 3   , N Colors
% CMYK:  N by 4
%

if isempty(c)
    return
end

s = size(c);
n = size(s,2);

m = s(2+(n==3));

sub = { ':'  ':' };
sub = sub(1:(1+(n==3)));

if m == 3
    %  RGB --> CMYK
    c = 1 - c;
    k = min(c,[],n);
    c = c - k(sub{:},[1 1 1]);
    c = cat(n ,c ,k );
else
    % CMYK --> RGB
    c = 1 - ( c(sub{:},[1 2 3]) + c(sub{:},[4 4 4]) );
end

end