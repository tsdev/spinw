function [TR, S] = orbmesh(qNum, varargin)
% creates triangular mesh of selected hydrogen orbital isosurfaces
%
% [TR,S] = SWPLOT.ORBMESH(qNum, 'Option1', Value1, ...)
%
% [TR,S] = SWPLOT.ORBMESH(qLabel, 'Option1', Value1, ...)
%
% The output is a triangulated surface of the selected orbital, containing
% ~nPatch^3 triangular faces. The output can be plotted using the trimesh()
% function. Using the qLabel input, all orbitals corresponds to the same
% probability isosurface (P ~ 1e-5), thus the relative sizes of the s, p
% and d orbitals are correct.
%
% Input:
%
% qNum      Quantum numbers in vector (n, l, m, {pm}), where pm defines the
%           optional linear combination of the +m and -m orbitals:
%               PSI = PSI(n,l,m) + pm*PSI(n,l,-m)
%           The default value of pm is 0. If pm is +1 or -1, real orbitals
%           will be produced.
% qLabel    Label of the orbital in a string, corresponding quantum
%           numbers:
%               's'         qNum = [1 0 0  0]
%               'p_x'       qNum = [2 1 1  1]
%               'p_y'       qNum = [2 1 0 -1]
%               'p_z'       qNum = [2 1 0  0]
%               'd_xy'      qNum = [3 2 2 -1]
%               'd_xz'      qNum = [3 2 1  1]
%               'd_yz'      qNum = [3 2 1 -1]
%               'd_z2'      qNum = [3 2 0  0]
%               'd_x2-y2'   qNum = [3 2 2  1]
%
% Options:
%
% rLim      Limit of the axes, default is 8 Angstrom. It determines the
%           probability value where the isosurface is drawn.
% nPatch    Number of points in the surface mesh along the three
%           dimensions, default is 20.
%
% Output:
%
% TR        TriRep class triangulation object for plotting with trimesh().
% S         Color index per face, the value is 1 or 2 for negative and
%           positive wave function values. Can be used to differentiate
%           between the negative and positive wave vector values with
%           colors.
%

if nargin == 0
    help swplot.orbmesh
    return
end

% Bohr radius
rBohr = 0.52917721067; % Angstrom
inpForm.fname  = {'nPatch' 'rLim' 'P'  };
inpForm.defval = {20       []     []  };
inpForm.size   = {[1 1]    [1 1]  [1 1]};
inpForm.soft   = {false    true   true };

param = sw_readparam(inpForm,varargin{:});

if ischar(qNum)
    switch qNum
        case 'd_xy'
            qNum = [3 2 2 -1];
            rLim = 8.0;
        case 'd_xz'
            qNum = [3 2 1 1];
            rLim = 8.0;
        case 'd_yz'
            qNum = [3 2 1 -1];
            rLim = 8.0;
        case 'd_z2'
            qNum = [3 2 0  0];
            rLim = 8.1;
        case 'd_x2-y2'
            qNum = [3 2 2  1];
            rLim = 9.5;
        case  'p_x'
            qNum = [2 1 1  1];
            rLim = 6.0;
        case 'p_y'
            qNum = [2 1 1 -1];
            rLim = 6.0;
        case 'p_z'
            qNum = [2 1 0  0];
            rLim = 5.0;
        case 's'
            qNum = [1 0 0  0];
            rLim = 2.3;
        otherwise
            error('orbmesh:WrongQNum','Wrong orbital label!');
    end
    
    if isempty(param.rLim)
        param.rLim = rLim;
    end
    
elseif isempty(param.rLim)
    param.rLim = 8;
end

% quantum numbers
n = qNum(1);
l = qNum(2); %  0<= l <n
m = qNum(3); % -l<= m <=l

if n<1
    error('orbmesh:WrongN','Wrong n quantum number (n>=1)!');
end
if (l<0) || (l>=n)
    error('orbmesh:WrongL','Wrong l quantum number (l>=0 & l<n)!');
end
if (m<-l) || (m>l)
    error('orbmesh:WrongM','Wrong m quantum number (-l<=m<=l)!');
end

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
        error('orbmesh:WrongInput','The pm value should be {-1,0,1}.')
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

% Sign the orbital,  - and + wave function sign
S = sign(psi);

% Electron density
psi = abs(psi).^2;

if isempty(param.P)
    % Find the probability iso value that touches the given plot range
    psi0 = psi;
    psi0(2:(end-1),2:(end-1),2:(end-1)) = 0;
    P = max(psi0(:));
else
    P = param.P;
end

[F,V,S] = isosurface(x,y,z,psi,P,S);

% create color selector
S = (S+3)/2;

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