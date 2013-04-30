function [w, s] = horace(obj, qh, qk, ql, p, varargin) %#ok<INUSL>
% dispersion/correltion function calculator, can be called from Horace
%
% [w, s] = HORACE(obj, qh, qk, ql, p, 'Option1, Value1, ...) function to produce
% spin wave dispersion and intensity for Horace (<a href=http://horace.isis.rl.ac.uk>http://horace.isis.rl.ac.uk</a>).
%
% Input:
%
% ojb           Input sw object.
% qh, qk, ql    Reciprocal lattice components in reciprocal lattice units.
% p             Prameters, not used.
%
% Options:
%
% swfunc        Which function to evaluate for spin wave calculation.
%               Default is @spinwave.
%
% Example:
%
% horace_on;
% d3dobj = d3d(cryst.abc,[0 1 0 0],[0,0.01,1],[0 0 1 0],[0,0.01,1],[0 0 0 1],[0,0.1,10]);
% d3dobj = disp2sqw_eval(d3dobj,cryst.horace,[],0.1);
% plot(d3dobj);
%
% This example creates a d3d object, a square in (h,k,0) plane and in
% energy between 0 and 10 meV. Then calculates the neutron scattering
% intensity of the spin wave in this volume and plots it using sliceomatic.
%
% See also SW, SW.SWINC, SW.SPINWAVE.
%

if nargin <= 1
    help sw.horace;
    return;
end

inpForm.fname  = {'swfunc'  };
inpForm.defval = {@spinwave };
inpForm.size   = {[1 1]     };

param   = sw_readparam(inpForm,varargin{:});

spectra = param.swfunc(obj,[qh(:) qk(:) ql(:)]',param);
spectra = sw_neutron(spectra,'pol',false);

nModes = size(spectra.omega,1);
nHkl   = size(spectra.omega,2);

% dispersion in cell
w     = mat2cell(spectra.omega',nHkl,ones(nModes,1));

% intensity in cell
s     = mat2cell(spectra.Sperp' ,nHkl,ones(nModes,1));

end