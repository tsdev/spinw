function C = sw_bose(oldT,newT,E)
% coefficient for boson correlation functions for different temperatures
%
% C = SW_BOSE(oldT,newT,E)
%
% Input:
%
% oldT      Original temperature in Kelvin.
% newT      New temperature in Kelvin.
% E         Energy in meV, positive is the particle creation side (neutron
%           energy loss side in scattering experiment).
%
% Output:
%
% C         Correction coefficients that multiplies the correlation
%           function. If any of the input is a vector, C will be also a
%           vector with the same dimensions.
%

if nargin == 0
    help sw_bose
    return
end

kB   = 0.086173324; %    Boltzmann constant: k_B [meV/K]

% Bose factor for the original temperature
oldN = 1./(exp(abs(E)/(kB*oldT))-1) + double(E>0);
oldN(isinf(oldN)) = 1;
% Bose factor for the new temperature
newN = 1./(exp(abs(E)/(kB*newT))-1) + double(E>0);
newN(isinf(newN)) = 1;
% conversion factor
C = newN./oldN;

end
