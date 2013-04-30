function hSurf = sw_orbital(qNum, varargin)
% hSurf = SW_ORBITAL(qNum, 'Option1', Value1, ...) returns a polygon of the
% hydrogen orbitals.
%
% Input:
%
% qNum      Quantum numbers in vector (n, l, m, real). If the length of the
%           vector is 4, real orbitals will be produced, if real==1, the
%           sum of the +m and -m orbitals is calculated, however if
%           real==2, the difference is calculated.
%
% Options:
%
% surfRes   Number of points in the surface mesh along every axis,
%           default is 30.
% rLim      Limits of the axes, default is 32.
% P         Constant probability surface, default is 1E-5.
% rBohr     Bohr radius, default is 1.
% norm      Whether to normalise the axes, default is true.
%
% Radius is normalised to 1.
%
%

inpForm.fname  = {'surfRes' 'rLim' 'P'   'rBohr' 'norm'};
inpForm.defval = {30        32     1e-5  1       true  };
inpForm.size   = {[1 1]     [1 1]  [1 1] [1 1]   [1 1] };

param = sw_readparam(inpForm,varargin{:});


% quantum numbers
n = qNum(1);
l = qNum(2); %  0<= l <n
m = qNum(3); % -l<= m <=l

if n<1
    error('sw:sw_orbital:WrongN','Wrong n quantum number!');
end

if (l<0) || (l>=n)
    error('sw:sw_orbital:WrongL','Wrong l quantum number!');
end
if (m<-l) || (m>l)
    error('sw:sw_orbital:WrongM','Wrong m quantum number!');
end

% Bohr radius
rBohr = param.rBohr;
% normalization
N = abs(sign(m)*sqrt(2)+(sign(abs(m))-1)*2);
% angular part
SphericalYlm = @(l,m,theta,phi) sqrt((2*l+1)/(4*pi)*factorial(l-abs(m))/...
    factorial(l+abs(m)))*AssociatedLegendre(l,m,cos(theta)).*exp(1i*m*phi);

if length(qNum) == 4
    
    if qNum(4) == 1
        oSign = 1;
    else
        oSign = 1i;
    end
    
    Y = @(l,m,theta,phi) (SphericalYlm(l,m,theta,phi)+oSign^2*SphericalYlm(l,-m,theta,phi))/N/oSign;
else
    Y = @(l,m,theta,phi) SphericalYlm(l,m,theta,phi)/N;
end
% radial part
R = @(n,l,r) sqrt((2/(rBohr*n))^3*factorial(n-l-1)/(2*n*factorial(n+l))).*...
    exp(-r/(rBohr*n)).*(2*r/(rBohr*n)).^l*1/factorial(n-l-1+2*l+1).*...
    AssociatedLaguerre(n-l-1,2*l+1,2*r/(rBohr*n));
% wave function
psi = @(n,l,m,r,theta,phi) R(n,l,r).*Y(l,m,theta,phi);

% Setting the grid
rLim  = param.rLim;
surfRes = param.surfRes;
[x,y,z] = ndgrid(linspace(-rLim,rLim,surfRes),linspace(-rLim,rLim,surfRes),linspace(-rLim,rLim,surfRes));

% conversion Cartesian to spherical coordinates
r     = sqrt(x.^2+y.^2+z.^2);
theta = acos(z./r);
phi   = atan2(y,x);

% plot orbital,  - and + wave function phase
colors = sign(psi(n,l,m,r,theta,phi));

hSurf = isosurface(x,y,z,psi(n,l,m,r,theta,phi).^2,param.P,colors);

if param.norm
    rNorm = max(max(max(abs(hSurf.vertices))));
    hSurf.vertices = hSurf.vertices/rNorm;
end

end


function[Anm]=AssociatedLaguerre(n,m,x)
% Associated Laguerre

Anm = 0;
for ii = 0:n
    Anm = Anm+factorial(m+n)*nchoosek(m+n,n-ii)/factorial(ii)*(-x).^ii;
end
end

function[Alm]=AssociatedLegendre(l,m,x)
% Associated Legendre

Alm = 0;
for r = 0:floor(1/2*l-1/2*abs(m))
    Alm = Alm+(-1)^r*nchoosek(l-2*r,abs(m))*nchoosek(l,r)*nchoosek(2*l-2*r,l)*x.^(l-2*r-abs(m));
end
Alm = (-1)^m*(1-x.^2).^(abs(m)/2).*(factorial(abs(m))/2^l*Alm);
end