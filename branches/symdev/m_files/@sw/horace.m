function [w, s] = horace(obj, qh, qk, ql, ~)
% dispersion/correltion function calculator, can be called from Horace
%
% [w, s] = HORACE(obj, qh, qk, ql, p) function to produce
% spin wave dispersion and intensity for Horace (<a href=http://horace.isis.rl.ac.uk>http://horace.isis.rl.ac.uk</a>).
%
% Input:
%
% ojb           Input sw object.
% qh, qk, ql    Reciprocal lattice components in reciprocal lattice units.
% p             Prameters, not used.
%
% Example:
%
% horace_on;
% d3dobj = d3d(cryst.abc,[0 1 0 0],[0,0.01,1],[0 0 1 0],[0,0.01,1],[0 0 0 1],[0,0.1,10]);
% d3dobj = disp2sqw_eval(d3dobj,@cryst.horace,[],0.1);
% plot(d3dobj);
%
% This example creates a d3d object, a square in (h,k,0) plane and in
% energy between 0 and 10 meV. Then calculates the neutron scattering
% intensity of the spin wave in this volume and plots it using sliceomatic.
%
% See also SW, SW.SPINWAVE.
%

if nargin <= 1
    help sw.horace;
    return;
end

<<<<<<< .mine
if nargin > 5
    spectra = obj.spinwave([qh(:) qk(:) ql(:)]',varargin{:});
else
    spectra = obj.spinwave([qh(:) qk(:) ql(:)]');
end

=======
spectra = obj.spinwave([qh(:) qk(:) ql(:)]');
>>>>>>> .r62
spectra = sw_neutron(spectra,'pol',false);

% add all modes for different twins
if iscell(spectra.omega)
    omega = cell2mat(spectra.omega');
    Sperp = cell2mat(spectra.Sperp');
else
    omega = spectra.omega;
    Sperp = spectra.Sperp;
end

nModes = size(omega,1);
nHkl   = size(omega,2);

% dispersion in cell
w = mat2cell(omega',nHkl,ones(nModes,1));
% intensity in cell
s = mat2cell(Sperp' ,nHkl,ones(nModes,1));

end