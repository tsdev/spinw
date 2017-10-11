function out = sw_converter(value, unitIn, unitOut, particleName,invert)
% converts energy and momentum units for a given particle
% 
% ### Syntax
% 
% `out = sw_converter(value, unitIn, unitOut)`
% 
% `out = sw_converter(value, unitIn, unitOut, particleName)`
%
% ### Description
% 
% `out = sw_converter(value, unitin, unitout)` will convert momentum and
% energy values assuming neutron as a particle.
% 
% `out = sw_converter(value, unitin, unitout,particleName)` will convert
% momentum and energy values for a given particle, such as neutron, photon,
% etc.
%
% ### Example
%
% Calculate the energy of a neutron (in meV) which has a wavelength of
% 5 \\Angstrom:
%
% ```
% >>sw_converter(5,'A','meV')>>
% ```
%
% Calculate the wavelength of X-ray in \\Angstrom that has 7.5 keV energy:
%
% ```
% >>sw_converter(7.5,'keV','A','photon')>>
% ```
%
% ### Input Arguments
% 
% `value`
% : Numerical input value, can be scalar or matrix with arbitrary
%   dimensions.
%
% `unitIn`
% : Units of the input value, one of the following string:
%   * `'A-1'`        momentum in \\Angstrom$^{-1}$,
%   * `'k'`          momentum in \\Angstrom$^{-1}$,
%   * `'Angstrom'`   wavelength in \\Angstrom,
%   * `'lambda'`     wavelength in \\Angstrom,
%   * `'A'`          wavelength in \\Angstrom,
%   * `'K'`          temperature in Kelvin,
%   * `'m/s'`        speed in m/s,
%   * `'J'`          energy in Joule,
%   * `'meV'`        energy in meV,
%   * `'eV'`         energy in eV,
%   * `'keV'`        energy in keV,
%   * `'THz'`        frequency in Thz,
%   * `'cm-1'`       $2\pi/\lambda$ in cm$^{-1}$,
%   * `'fs'`         wave period time in fs,
%   * `'ps'`         wave period time in ps,
%   * `'nm'`         wavelength in nm,
%   * `'um'`         wavelength in $\mu$m.
% 
% `unitOut`
% : Units of the output value, same strings are accepted as for `unitIn`.
%
% `particleName`
% : String, the name of the particle, one of the following values: 
%   `'neutron'` (default), `'proton'`, `'electron'`, `'photon'`, `'xray'`,
%   `'light'`.
%

if nargin == 0
    help sw_converter
    return
end

if isempty(unitIn) && isempty(unitOut)
    out = value;
    return
elseif isempty(unitIn) || isempty(unitOut)
    error('sw_converter:WrongInput','One of the unit is empty string!')
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
if nargin < 4
    particleName = 'neutron';
end

if nargin < 5
    invert = false;
end

% flip the units
if invert
    unitTemp = unitIn;
    unitIn   = unitOut;
    unitOut  = unitTemp;
end

switch particleName
    case 'neutron'
        m = mN;
    case 'proton'
        m = mP;
    case 'electron'
        m = mE;
    case {'photon' 'xray' 'x-ray' 'light'}
        m = 0;
    otherwise
        error('sw_converter:WrongParticle','Particle name is wrong!');
end

% Conversions
switch unitIn
    % convert everything into momentum in Angstrom^-1
    case {'Angstrom' 'A' 'lambda' 'angstrom'}
        k = 2*pi./value;
    case 'nm'
        k = 2*pi./(value*10);
    case 'um'
        k = 2*pi./(value*1e4);
    case {'Kelvin' 'K'}
        if m~=0
            k = sqrt((value*EK2J*2*m))/hBar/1e10;
        else
            k = value*EK2J/clight/1e10/hBar;
        end
    case {'k' 'A-1' 'Angstrom-1' 'angstrom-1'}
        k = value;
    case 'mps'
        if m~=0
            k = value*m/hBar/1e10;
        else
            error('sw_converter:WrongUnit','Speed cannot be an input for photon!');
        end
    case 'meV'
        if m~=0
            k = sqrt((value*EmeV2J*2*m))/hBar/1e10;
        else
            k = value/hBar/clight*EmeV2J/1e10;
        end
    case 'J'
        if m~=0
            k = sqrt((value*2*m))/hBar/1e10;
        else
            k = value/hBar/clight/1e10;
        end
        
    case 'eV'
        if m~=0
            k = sqrt((value*1e3*EmeV2J*2*m))/hBar/1e10;
        else
            k = value*1e3/hBar/clight*EmeV2J/1e10;
        end
    case 'keV'
        if m~=0
            k = sqrt((value*1e6*EmeV2J*2*m))/hBar/1e10;
        else
            k = value*1e6/hBar/clight*EmeV2J/1e10;
        end
    case 'THz'
        if m~=0
            k = sqrt(value*ETHz2meV*EmeV2J*2*m)/hBar/1e10;
        else
            k = value*1e12*2*pi/clight/1e10;
        end
    case 'cm-1'
        k = 1e-8*value*2*pi;
    case 'fs'
        k = 1/(value*1e-15)*2*pi/clight/1e10;
    case 'ps'
        k = 1/(value*1e3*1e-15)*2*pi/clight/1e10;
end

switch unitOut
    % convert from momentum in Angstrom^-1 to output units
    case {'Angstrom' 'A' 'angstrom' 'lambda'}
        out = 2*pi./k;
    case 'nm'
        out = 2*pi./k*0.1;
    case 'um'
        out = 2*pi./k*1e-4;
    case {'Kelvin' 'K'}
        if m~=0
            out = (k*hBar*1e10).^2/2/m/EK2J;
        else
            out = k/EK2J*clight*1e10*hBar;
        end
    case {'k' 'A-1' 'Angstrom-1' 'angstrom-1'}
        out = k;
    case {'mps' 'm/s'}
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
    case 'eV'
        if m~=0
            out = (k*hBar*1e10).^2/EmeV2J/2/m*1e-3;
        else
            out = k/EmeV2J*clight*1e10*hBar*1e-3;
        end
    case 'keV'
        if m~=0
            out = (k*hBar*1e10).^2/EmeV2J/2/m*1e-6;
        else
            out = k/EmeV2J*clight*1e10*hBar*1e-6;
        end
    case 'THz'
        if m~=0
            out = (k*hBar*1e10).^2/ETHz2meV/EmeV2J/2/m;
        else
            out = k*clight/2/pi/1e12*1e10;
        end
    case 'cm-1'
        out = 1e8*k/2/pi;
    case 'fs'
        out = 2*pi/clight/k/1e10*1e15;
    case 'ps'
        out = 2*pi/clight/k/1e10*1e15*1e-3;
    case {'J' 'Joule'}
        if m~=0
            %k = sqrt((value*2*m))/hBar/1e10;
            out = (k*1e10*hBar)^2/2/m;
        else
            %k = value/hBar/clight/1e10;
            out = k*hBar*clight*1e10;
        end
end

end