function out = sw_converter(unitIn, value, unitOut, particleName)
% out = SW_CONVERTER(unitIn, value, unitOut, {particleName}) converts energy
% and momentum units for a given particle.
%
% Input:
%
% particleName      Name of the particle:
%                       'neutron'   default
%                       'proton'
%                       'electron'
%                       'photon'
% unitIn            Units of the input value:
%                       'k'         momentum in Angstrom^-1.
%                       'Angstrom'  wavelength in Angstrom.
%                       'Kelvin'    temperature in Kelvin.
%                       'mps'       speed in m/s.
%                       'meV'       energy in meV.
%                       'THz'       energy in Thz.
% unitOut           Units of the output value, same options as for unitIn.
% value             Numerical input value, can be arbitrary matrix.
%

if nargin == 0
    help sw_converter;
    return;
end

% constants in SI units
% Boltzmann constant
kB           = 1.3806488e-23;
% Planck constant in J*s
h            = 6.626068e-34;
% h-bar in J*s
hBar         = h/2/pi;
% electron charge in Coulomb
e            = 1.60217646e-19;
% neutron mass in kg
mN           = 1.674927351e-27;
% proton mass
mP           = 1.672621777e-27;
% electron mass
mE           = 9.10938291e-31;
% speed of light in m/s
clight       = 299792458;

% conversion: meV --> J
EmeV2J = e/1000;
% conversion: K --> J
EK2J = kB;
% conversion: meV --> THz
EmeV2THz = EmeV2J/hBar*1e-12/2/pi;
% conversion: THz --> meV
ETHz2meV = 1/EmeV2THz;

% Choose particle my selecting its mass
if nargin == 3
    particleName = 'neutron';
end

switch particleName
    case 'neutron'
        m = mN;
    case 'proton'
        m = mP;
    case 'electron'
        m = mE;
    case 'photon'
        m = 0;
    otherwise
        error('sw:sw_converter:WrongParticle','Particle name is wrong!');
end

% Conversions
switch unitIn
    % convert everything into momentum in Angstrom^-1
    case 'Angstrom'
        k = 2*pi./value;
    case 'Kelvin'
        if m~=0
            k = sqrt((value*EK2J*2*m))/hBar/1e10;
        else
            k = value*EK2J/clight/1e10;
        end
    case 'k'
        k = value;
    case 'mps'
        if m~=0
            k = value*m/hBar/1e10;
        else
            error('sw:sw_converter:WrongUnit','Speed cannot be an input for photon!');
        end
    case 'meV'
        if m~=0
            k = sqrt((value*EmeV2J*2*m))/hBar/1e10;
        else
            k = value*EmeV2J/clight/1e10;
        end
    case 'THz'
        if m~=0
            k = sqrt(value*ETHz2meV*EmeV2J*2*m)/hBar/1e10;
        else
            k = value*h*1e12/1e10;
        end
end

switch unitOut
    % convert from momentum in Angstrom^-1 to output units
    case 'Angstrom'
        out = 2*pi./k;
    case 'Kelvin'
        if m~=0
            out = (k*hBar*1e10).^2/2/m/EK2J;
        else
            out = k/EK2J*clight*1e10;
        end
    case 'k'
        out = k;
    case 'mps'
        if m~=0
            out = k/m*hBar*1e10;
        else
            out = clight;
        end
    case 'meV'
        if m~=0
            out = (k*hBar*1e10).^2/EmeV2J/2/m;
        else
            out = k/EmeV2J*clight*1e10*hBar;
        end
    case 'THz'
        if m~=0
            out = (k*hBar*1e10).^2/ETHz2meV/EmeV2J/2/m;
        else
            out = k/h/1e12*1e10;
        end
end

end