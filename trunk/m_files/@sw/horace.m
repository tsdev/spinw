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
% See also SW, SW.SWINC, SW.SPINWAVE.
%

if nargin <= 1
    help sw.horace;
    return;
end

inpForm.fname  = {'swfunc'  };
inpForm.defval = {@spinwave };
inpForm.size   = {[1 1]     };

param  = sw_readparam(inpForm,varargin{:});

spec   = param.swfunc(obj,[qh(:) qk(:) ql(:)]',param);
spec   = sw_neutron(spec,'pol',false);

nModes = size(spec.omega,1);
nHkl   = size(spec.omega,2);

% dispersion in cell
w     = mat2cell(omega',nHkl,ones(nModes,1));

% intensity in cell
s     = mat2cell(spec.Sperp' ,nHkl,ones(nModes,1));

end