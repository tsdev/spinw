function Qm = sw_econtract(Q,varargin)
% converts (Q,E) values to Q values for diffraction instrument
% 
% ### Syntax
% 
% `Qm = sw_econtract(Q,Name,Value)`
% 
% ### Description
% 
% `Qm = sw_econtract(Q,Name,Value)` converts $(Q,E)$ values in the phase
% space into $Q$ values as it would appear when measuring it with via
% neutron diffraction.
% 
% ### Input Arguments
% 
% `Q`
% : Input values in reciprocal space in the scattering plane in
%   \\ang$^{-1}$ units, dimensions are $[2\times n_Q]$.
% 
% ### Name-Value Pair Arguments
% 
% `'E'`
% : Energy transfer value in meV, default value is zero.
% 
% `'lambda'`
% : Wavelength of the incident neutron beam in \\ang.
% 
% `'ki'`
% : Momentum of the incidend neutron beam in \\ang$^{-1}$, alternative
%   input to `lambda`.
% 
% `'sense'`
% : Scattering sense:
%
%   * `1`  detectors are on the right hand side from the incident beam direction, default.
%   * `-1` detectors are on the left hand side from the incident beam direction.
% 
% ### See Also
% 
% [sw_converter]
%

if nargin == 0
    swhelp sw_econtract
    return
end

inpForm.fname  = {'E'   'lambda' 'ki'  'sense'};
inpForm.defval = {0     0        0     1      };
inpForm.size   = {[1 1] [1 1]    [1 1] [1 1]  };

param = sw_readparam(inpForm,varargin{:});

if param.lambda == 0
    % Angstrom
    ki = param.ki;
else
    % Angstrom^-1
    ki = 2*pi/param.lambda;
end

if ki == 0
    error('sw_flatcone:kierror','Missing ki or lambda!');
end

% Q absolute values
Qabs = sqrt(sum(Q.^2,1));
% conversion from meV to Angstrom^-1 for neutrons
meV2k = sw_converter('meV',1,'k','neutron');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (Q, omega) --> (kf, theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% energy loss due to magnons
kf = sqrt(ki^2-param.E*meV2k);
% theta scattering angle
theta = acos((ki^2 + kf^2 - Qabs.^2)/2/ki/kf)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (kf, theta) --> (Qmeas, phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% new momentum transfer absolute value
Qmabs = 2*ki*sin(theta);
% rotation of the momentum
phi = acos((Qabs.^2 + Qmabs.^2 - (ki-kf).^2)./(Qabs.*Qmabs)/2);
% remove small imaginary component
phi(abs(imag(phi))<1e-6) = real(phi(abs(imag(phi))<1e-6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate Qmeas vector in Angstrom^-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% phi value of the original momentum transfer
phi0 = atan2(Q(2,:),Q(1,:));
% new momentum transfer vector
phim = phi0+param.sense*phi;
Qm = bsxfun(@times,Qmabs,[cos(phim); sin(phim)]);

% select where the triangle cannot be closed
Qm(:,isnan(sum(Qm,1))) = 0;
Qm(:,logical(imag(sum(Qm,1)))) = 0;
Qm = real(Qm);

end