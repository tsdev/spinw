function [TR,C] = sw_orbital(qNum, varargin)
% creates triangular mesh of selected hydrogen orbital isosurfaces
%
% [TR,C] = ORBITAL(qNum, 'Option1', Value1, ...)
%
% The output is a triangulated surface of the selected orbital, containing
% ~nPatch^3 triangular faces. The output can be plotted using the trimesh()
% function.
%
% Input:
%
% qNum      Quantum numbers in vector (n, l, m, {pm}), where pm defines the
%           optional linear combination of the +m and -m orbitals:
%               PSI = PSI(n,l,m) + pm*PSI(n,l,-m)
%           The default value of pm is 0. If pm is +1 or -1, real orbitals
%           will be produced.
%
% Options:
%
% nPatch    Number of points in the surface mesh along the three
%           dimensions, default is 20.
% rLim      Limits of the axes, default is 20, determines the probability
%           value where the isosurface is drawn. Default value is 20.
% rBohr     Bohr radius, default is 1.
% norm      Whether to normalise the axes to 1, default is true.
% color     Colors for positive and negative wave function values. Either a
%           cell of two strings (color names, see swplot.color) or 2 RGB
%           color values in a 3x2 matrix. Default is {'blue' 'red'}.
%
% Output:
%
% TR        TriRep class triangulation object for plotting with trimesh().
% C         Color data for each triangular faces (blue for 

if nargin == 0
    help sw_orbital
    return
end

C0 = [0 0 1;1 0 0]';

inpForm.fname  = {'nPatch' 'rLim' 'rBohr' 'norm' 'color'};
inpForm.defval = {20       20     1       true   C0     };
inpForm.size   = {[1 1]    [1 1]  [1 1]   [1 1]  [-1 2] };

param = sw_readparam(inpForm,varargin{:});

% quantum numbers
n = qNum(1);
l = qNum(2); %  0<= l <n
m = qNum(3); % -l<= m <=l

if n<1
    error('orbital:WrongN','Wrong n quantum number!');
end
if (l<0) || (l>=n)
    error('orbital:WrongL','Wrong l quantum number!');
end
if (m<-l) || (m>l)
    error('orbital:WrongM','Wrong m quantum number!');
end

% Bohr radius
rBohr = param.rBohr;
% normalization
N = abs(sign(m)*sqrt(2)+(sign(abs(m))-1)*2);
% angular part
SphericalYlm = @(l,m,theta,phi) sqrt((2*l+1)/(4*pi)*factorial(l-abs(m))/...
    factorial(l+abs(m)))*AssociatedLegendre(l,m,cos(theta)).*exp(1i*m*phi);

if numel(qNum) == 3
    qNum(4) = 0;
end

% linear combination
pm = qNum(4);

if pm == 0
    % imaginary orbitals for m~=0
    Y = @(l,m,theta,phi) SphericalYlm(l,m,theta,phi)/N;
else
    % real orbitals
    if pm == 1
        oSign = 1;
    elseif pm == -1
        oSign = 1i;
    else
        error('orbital:WrongInput','The pm value should be {-1,0,1}.')
    end
    
    Y = @(l,m,theta,phi) (SphericalYlm(l,m,theta,phi)+oSign^2*SphericalYlm(l,-m,theta,phi))/N/oSign;
end
% radial part
R = @(n,l,r) sqrt((2/(rBohr*n))^3*factorial(n-l-1)/(2*n*factorial(n+l))).*...
    exp(-r/(rBohr*n)).*(2*r/(rBohr*n)).^l*1/factorial(n-l-1+2*l+1).*...
    AssociatedLaguerre(n-l-1,2*l+1,2*r/(rBohr*n));
% wave function
psi = @(n,l,m,r,theta,phi) R(n,l,r).*Y(l,m,theta,phi);

% Setting the grid
rl      = param.rLim;
nP      = param.nPatch;
rV      = linspace(-rl,rl,nP);
[x,y,z] = ndgrid(rV,rV,rV);

% conversion Cartesian to spherical coordinates
r     = sqrt(x.^2+y.^2+z.^2);
theta = acos(z./r);
phi   = atan2(y,x);
% Calculate the function on the grid
psi   = psi(n,l,m,r,theta,phi);

% Color the orbital,  - and + wave function sign
sCol = sign(psi);

% Electron density
psi = psi.^2;

% Find the probability iso value that touches the given plot range
psi0 = psi;
psi0(2:(end-1),2:(end-1),2:(end-1)) = 0;
P = max(psi0(:));

[F,V,C] = isosurface(x,y,z,psi,P,sCol);

if param.norm
    V = V/max(V(:));
end

% Face color data
col = swplot.color(param.color)';
C   = col((C+3)/2,:)/255;

% create triangulation
TR  = TriRep(F,V); %#ok<DTRIREP>

end


function Anm = AssociatedLaguerre(n,m,x)
% Associated Laguerre

Anm = 0;
for ii = 0:n
    Anm = Anm+factorial(m+n)*nchoosek(m+n,n-ii)/factorial(ii)*(-x).^ii;
end

end

function Alm = AssociatedLegendre(l,m,x)
% Associated Legendre

Alm = 0;
for r = 0:floor(1/2*l-1/2*abs(m))
    Alm = Alm+(-1)^r*nchoosek(l-2*r,abs(m))*nchoosek(l,r)*nchoosek(2*l-2*r,l)*x.^(l-2*r-abs(m));
end
Alm = (-1)^m*(1-x.^2).^(abs(m)/2).*(factorial(abs(m))/2^l*Alm);

end