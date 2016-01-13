function weight = horace_sqw(obj, qh, qk, ql, en, p_in)
% calculates spin wave dispersion/correlation functions to be called from Horace
%
% weight = HORACE_SQW(obj, qh, qk, ql, en, p)
%
% The function produces spin wave dispersion and intensity for Horace
% (<a href=http://horace.isis.rl.ac.uk>http://horace.isis.rl.ac.uk</a>).
%
% Input:
%
% obj        Input sw object.
% qh, qk, ql Reciprocal lattice components in reciprocal lattice units.
% en         Energy transfers at which to calculate S(q,w)
% p          Parameters, in the order defined by horace_setpar.
%            In addition, if numel(p)>numel(mapping) where mapping is
%            the cell array defined in horace_setpar, the next 3 values
%            are taken to be:
%              I0   - amplitude to scale calculated S(q,w) by
%              fwhh - full-width at half height to convolve with a Gaussian
%              bkgd - constant background to add to calculated S(q,w)
%            Default values are: I0=1, fwhh=max(energy)/100, bkgd=0;
%
% Output:
%
% weight     Array of neutron intensity at specified hkl, energy points.
%
% Example:
%
% ...
% horace_on;
% tri = sw_model('triAF',2);
% tri.horace_setpar('mapping',{'J1' 'J2'},'fwhm',0.2);
% d3dobj = d3d(tri.abc,[0 1 0 0],[0,0.01,1],[0 0 1 0],[0,0.01,1],[0 0 0 1],[0,0.1,10]);
% d3dobj = sqw_eval(d3dobj,@cryst.horace_sqw,[1 0.5 1.5 0.5 0.01]);
% plot(d3dobj);
%
% This example creates a d3d object, a square in (h,k,0) plane and in
% energy between 0 and 10 meV. Then calculates the inelastice neutron
% scattering intensity of triagular lattice antiferromagnet and plots it
% using sliceomatic.
%
% See also SW, SW.SPINWAVE, SW.MATPARSER, SW.HORACE_SETPAR, SW_READPARAM.
%

if nargin <= 1
    help sw.horace_sqw;
    return;
end

nPar = numel(obj.matrix.horace.mapping);
I0 = 1;
bkgd = 0;
if numel(p_in) > nPar
    I0 = p_in(nPar+1);
end
if numel(p_in) > (nPar+1)
    fwhh = p_in(nPar+2);
end
if numel(p_in) > (nPar+2)
    bkgd = p_in(nPar+3);
end

[e, sf] = obj.horace(qh,qk,ql,p_in(1:nPar));

if ~exist('fwhh','var') || ~isnumeric(fwhh)
    fwhh = max(cellfun(@max,w))/100;
end

sig = fwhh/sqrt(log(256));
weight = zeros(numel(qh),1);
for ii=1:numel(e)
    weight = weight + sf{ii}(:).*exp(-(e{ii}(:)-en(:)).^2/(2*sig^2))/(sig*sqrt(2*pi));
end    
weight = I0*reshape(weight,size(qh)) + bkgd;

end
