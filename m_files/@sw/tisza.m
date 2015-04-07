function optRes = tisza(obj,varargin)
% find magnetic structure using Luttinger-Tisza method
%
% optRes = TISZA(obj,'option1', value1 ...)
%
% The function determines the magnetic structure using the Luttinger-Tisza
% method.
%
% Options:
%
% n         Normal vector of the spiral plane, default is [0 0 1].
% hkl       Points in momentum space for minimizing the energy.
%
% See also SW.OPTMAGSTR, SW.OPTMAGSTEEP.
%

inpForm.fname  = {'n'    };
inpForm.defval = {[0 0 1]};
inpForm.size   = {[1 3]  };

param = sw_readparam(inpForm, varargin{:});

% Fourier transform of the magnetic interactions
ft = fourier(obj.hkl);




end