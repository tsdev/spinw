function cMat = sw_multicolor(vMat, cMap, cLim, nCol, cflip)
% overlays monochrome maps into a single RGB map
% 
% ### Syntax
% 
% `cmat = sw_multicolor(vmat, cmap, clim, {ncol}, {flipud})`
% 
% ### Description
% 
% `cmat = sw_multicolor(vmat, cmap, clim, {ncol}, {flipud})` takes 2D
% matrices and overlays them and generating RGB color additively from the
% user defined colors correspond to each map. The function is used in
% [sw_plotspec] when multiple correlation functions are overlayed on the
% same plot.
% 
% ### Examples
% 
% In this example we create two intensity maps stored in the square
% matrices `A` and `B` (linearly changing intensity along $x$ and $y$ axes
% respectively, intensity ranging between -2 and 2). We plot these
% intensity maps by converting them to RGB colors using the inline function
% `rgbMap` and the Matlab built-in function `image`. We use `sw_multicolor`
% function to additively combine `A` and `B` and providing a color
% saturation value of 1 (and lowest value of -1). It is clearly visible on
% the resulting plot of `C` that it is white where both `A` and `B` has
% zero value (or below the lowest color value of -1) and it is red+green
% where both `A` and `B` are saturated.
%
% ```
% >>rgbMap  = @(mat,RGB,clim)bsxfun(@plus,ones(1,1,3),bsxfun(@times,(mat-clim(1))/diff(clim),permute(RGB(:)-1,[2 3 1])));
% >>
% >>red   = [1;0;0];
% >>green = [0;1;0];
% >>
% >>[A,B] = ndgrid(linspace(-2,2,501),linspace(-2,2,501));
% >>C = sw_multicolor(cat(3,A,B),[red green],[-1 1]);
% >>
% >>figure
% >>
% >>subplot(1,3,1)
% >>image(rgbMap(A,red,[-1 1]))
% >>title A
% >>
% >>subplot(1,3,2)
% >>image(rgbMap(B,green,[-1 1]))
% >>title B
% >>
% >>subplot(1,3,3)
% >>image(C)
% >>title C=A+B
% >>>pos = get(gcf,'Position')
% >>>pos(4) = round(pos(4)*0.4)
% >>>set(gcf,'Position',pos)
% >>snapnow
% ```
% 
% ### Input Arguments
% 
% `vMat`
% : Matrix that contains the input 2D intensity data, dimensions are
%   $[d_1\times d_2\times n_{plot}]$, where each intensity map has a
%   dimension of $[d_1\times d_2]$.
% 
% `cMap`
% : Defines the color map that maps intensity values within the `cLim`
%   limits to colors, can be the following types:
%   * `matrix`  Matrix of RGB colors, where each column
%               corresponds to an RGB triplet. The dimension of the matrix
%               is $[3\times n_{plot}]$ and the $i$th color corresponds to
%               the color of the $i$th intensity map in the `vMat` stack.
%   * `cell`    Cell of $n_{plot}$ colormap functions. For example
%               `{@copper @gray}`.
% 
% `cLim`
% : Defines the maximum and minimum intensity values that the given color
%   map will span. Values in vMat smaller than the minimum and larger than
%   the maximum will be shown with the minimum and maximum color in the
%   colormap respectively. `cLim` is a row vector with 2 elements..
% 
% `nCol`
% : Number of colors in the colormap. Optional, default value is
%   100.
% 
% `flipud`
% : If `true` the given colormaps are inverted. Optional, default value is
%   `false`.
% 
% ### Output Arguments
% 
% `cMat`
% : Matrix that contains the RGB image, with dimensions of $[d_1\times
%   d_2\times 3]$. The image can be shown using the `image` built-in Matlab
%   command.
%
% ### See Also
%
% [sw_plotspec]
%

if nargin < 3
    swhelp sw_multicolor
    return
end

if nargin == 3
    nCol = 100;
end

if nargin < 5
    cflip = false;
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
if iscell(cMap)
    for ii = 1:nPlot
        if cflip
            cMapGen(:,1:3,ii) = flipud(cMap{ii}(nCol));
        else
            cMapGen(:,1:3,ii) = cMap{ii}(nCol);
        end
    end
else
    for ii = 1:nPlot

        cMapGen(:,1:3,ii) = bsxfun(@plus,ones(1,3),bsxfun(@times,cMap(:,ii)'-1,linspace(0,1,nCol)'));
        
        if cflip
            cMapGen(:,1:3,ii) = flipud(cMapGen(:,1:3,ii));
        end
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