function cmap = cm_tasty(n, method)
% custom colormap with (some) tasty colors
%
% cmap = tasty({n},{method})
%
% Colormap made of black, violet, strawberry, pumpkin, royal yellow and
% white.
%
% Input:
%
% n         Number of colors, default is 121.
% method    Interpolation method, allowed values are:
%               'pchip'     Default, cubic interpolation.
%               'linear'    Linear interpolation between the given colors.
%
% Composed by Sandor Toth, 4.05.2016
%

% set the main colors in RGB format
RGB0 = [...
    0     0   0; ... % Black
    %76   40 130; ... % Violet (G&S)
    %    237  28  36; ... % Pigment red
    %194  63 255; ...
    155 50 204;...
    %    255 117  24; ... % Pumpkin
    %    250 218  94; ... % Royal yellow
    150 194 255; ...
    255 255 255];    % White
%    84   90 167; ... % Liberty
%    219 178 209; ... % Pink lavender
%    232  12  96; ... % stuff
%    210  14   7; ... % Strawberry
%    255 145 164; ... % Salmon (crayola)

RGB0 = RGB0/255;

% default colormap size 121x3
if nargin == 0
    n = 121;
end

if nargin < 2
    method = 'pchip';
end

% initial interpolation
n0 = 21;
% line integral
lInt = cumsum([0;sqrt(sum(diff(RGB0,1).^2,2))]);
% color positions along the line integral
cIdx = linspace(lInt(1),lInt(end),n0);

% cubic interpolation between the edge points
R = interp1(lInt,RGB0(:,1),cIdx,method);
G = interp1(lInt,RGB0(:,2),cIdx,method);
B = interp1(lInt,RGB0(:,3),cIdx,method);

RGB1 = [R;G;B]';

% convert to LAB color space where the color distance is well defined
LAB0 = rgb2lab(RGB1);

% line integral
lInt = cumsum([0;abs(diff(LAB0(:,1)))]);
%lInt = cumsum([0;sqrt(sum(diff(RGB0,1).^2,2))]);
% color positions along the line integral
cIdx = linspace(lInt(1),lInt(end),n);

% cubic interpolation between the edge points
R = interp1(lInt,RGB1(:,1),cIdx,method);
G = interp1(lInt,RGB1(:,2),cIdx,method);
B = interp1(lInt,RGB1(:,3),cIdx,method);

% convert back to RGB
cmap = [R;G;B]';

end