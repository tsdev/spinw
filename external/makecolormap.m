function cMap = makecolormap(varargin)
%% MAKECOLORMAP makes smoothly varying colormaps
% a = makeColorMap(beginColor, middleColor, endColor, numSteps);
% a = makeColorMap(beginColor, endColor, numSteps);
% a = makeColorMap(beginColor, middleColor, endColor);
% a = makeColorMap(beginColor, endColor);
%
% all colors are specified as RGB triples
% numSteps is a scalar saying howmany points are in the colormap
%
% Examples:
%
% peaks;
% a = makeColorMap([1 0 0],[1 1 1],[0 0 1],40);
% colormap(a)
% colorbar
%
% peaks;
% a = makeColorMap([1 0 0],[0 0 1],40);
% colormap(a)
% colorbar
%
% peaks;
% a = makeColorMap([1 0 0],[1 1 1],[0 0 1]);
% colormap(a)
% colorbar
%
% peaks;
% a = makeColorMap([1 0 0],[0 0 1]);
% colormap(a)
% colorbar
%
% Reference:
% A. Light & P.J. Bartlein, "The End of the Rainbow? Color Schemes for
% Improved Data Graphics," Eos,Vol. 85, No. 40, 5 October 2004.
% http://geography.uoregon.edu/datagraphics/EOS/Light&Bartlein_EOS2004.pdf
%

defaultNum = 128;

sz = cellfun('prodofsize',varargin);

constraits = vertcat(varargin{sz == 3});
steps = [varargin{sz == 1}];

numcstr = size(constraits,1);
numsteps = numel(steps);

if numsteps < 1
    steps = round(diff(linspace(1,defaultNum+1,numcstr)));
elseif numsteps == 1
    steps = round(diff(linspace(1,steps+1,numcstr)));
elseif numsteps > numcstr-1
    steps = steps(1:numcstr-1);
elseif numcstr-1 ~= numsteps
    steps = round(diff(linspace(1,defaultNum+1,numcstr)));
end

steps(1:end-1) = steps(1:end-1)+1;

cMap = [];
for k=1:numcstr-1
    cMap = [cMap(1:end-1,:); interpMap(constraits(k,:), constraits(k+1,:), steps(k))];
end

end

function cMap = interpMap(colorStart, colorEnd, n)
cMap = zeros(n,3);
for i = 1:3
    cMap(:,i) = linspace(colorStart(i), colorEnd(i), n);
end
end