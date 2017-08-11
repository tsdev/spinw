function setunit(obj, varargin)
% sets the physical units
%
% SETUNIT(obj, 'option1', value1 ...)
%
% Input:
%
% obj       SpinW object.
%
% Options:
%
% mode      Type of unit system, defined by one of the following strings:
%               'AmeVTK'    Typical units used in neutron/xray scattering:
%                               [Angstrom, meV, Tesla and Kelvin]
%               '1'         No units, all conversion factors are set to 1.
%

inpForm.fname  = {'mode'  };
inpForm.defval = {'AmeVTK'};
inpForm.size   = {[1 -1]  };

param = sw_readparam(inpForm, varargin{:});

switch param.mode
    case 'AmeVTK'
        k_B  = 0.086173324;     % Boltzmann constant: k_B [meV/K]
        mu_B = 0.057883818066;  % Bohr magneton: mu_B [meV/T]
        e    = 1.602176565e-19; % electron charge: e [C]
        mu0  = 4*pi*10*e;       % vacuum permeability: mu0 [T^2*A^3/meV]
        
        obj.unit.kB    = k_B;
        obj.unit.muB   = mu_B;
        obj.unit.mu0   = mu0;
        obj.unit.label = {char(197) 'meV' 'T' 'K'};
    case '1'
        obj.unit.kB    = 1;
        obj.unit.muB   = 1;
        obj.unit.mu0   = 1;
        obj.unit.label = {'' '' '' ''};
    otherwise
        error('spinw:setunit:WrongInput','The given ''mode'' option is invalid!')
end

end