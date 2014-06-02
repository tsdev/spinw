function spectra = sw_spinmotion(spectra, varargin)
% SW_SPINMOTION(spectra, 'option1', value1 ...) calculates the amplitude
% and phase of the local spin precessions belonging to a selected normal
% magnon mode.
%
% Options:
%
% iMagnon   Index of the normal magnon mode. Number has to be between one
%           and nMagExt, the number of spins in the magnetic supercell.
%           Default is the first magnon mode.
% Q         Momentum of the magnon in reciprocal lattice units (r.l.u.).
%           Either a vector (Qh,Qk,Ql) where spectra is calculated with
%           dimensions of [3 1], or scalar denoting the index of the Q
%           point in the list of calculated Q points. Default is the first
%           Q point.
%
% Output:
%
% Spectra contains an additional spectra.motion field with the following
% subfields:
%
% Amp1, Amp2    Orthogonal spin wave amplitude vectors for every spin in
%               the magnetic supercell, dimensions are [3 nMagExt]. The
%               selected normal magnon will create on the i-th spin a
%               precession around the ground state direction along the
%               ellipse defined by the Amp1(:,i) and Amp2(:,i) orthogonal
%               vectors as major and minor axes. At t=0 the spin deviation
%               vector points along Amp1.
% iQ            Index the magnon momentum, spectra.hkl(:,iQ) gives the
%               magnon momentum in r.l.u. units.
% iMagnon       The index of the normal magnon mode.
% omega         Energy of the selected normal magnon mode. It gives also
%
% See also SW.SPINWAVE.
%

inpForm.fname  = {'iMagnon' 'Q'   };
inpForm.defval = {1         1     };
inpForm.size   = {[1 1]    [-1 1] };

param = sw_readparam(inpForm, varargin{:});

if ~isfield(spectra,'T')
    error('sw_spinmotion:MissingField',['Input spectra has to be calculated '...
        'before using sw.spinwave() function with the ''saveT'' option set to true.']);
end

% sw object of the calculated spectrum
obj = spectra.obj;

iMagnon = param.iMagnon;

if numel(param.Q) == 1
    iQ = min(param.Q,size(spectra.hkl,2));
else
    if size(param.Q,1) ~= 3
        error('sw_spinmotion:WrongInput','Option ''Q'' has to be either column vector with 3 elements or a scalar!');
    end
    dQ = sqrt(sum((obj.basisvector*bsxfun(@minus,spectra.hkl,param.Q)).^2,1));
    iQ = find(dQ == min(dQ));
    
    if dQ(iQ) > 1e-3
        warning('sw_spinmotion:WrongQ','The closest calculated Q point (%5.3f %5.3f %5.3f) is far!',spectra.hkl(1,iQ),spectra.hkl(2,iQ),spectra.hkl(3,iQ));
    end
end

% number of magnetic atoms in the supercell
nMagExt = obj.nmagext;

if (iMagnon<1) || (iMagnon > nMagExt)
    error('sw_spinmotion:WrongInput',['The index of the normal magnon mode has to'...
        ' be between one and the number of moments in the magnetic supercell!']);
end

% select the positive dispersion modes
if spectra.omega(iMagnon,iQ)>0
    nMagS = iMagnon;
else
    nMagS = iMagnon + nMagExt;
end

% vector transforms the selected normal magnon mode into the localized
% magnons on the spins
Tm = spectra.T(:,nMagS,iQ);

% energy of the selected magnon
om = spectra.omega(nMagS,iQ);

% original complex prefactors of the boson operators
A = Tm(1:nMagExt);
B = Tm((nMagExt+1):end);

% complex amplitudes along the e1 and e2 axes (x,y)
Cx = transpose(A + B);
Cy = transpose(B - A);

% local spin coordinate system (e1, e2, e3||S)
% spin precession is in the (e2,e3) plane
M0 = obj.magtable;
S0 = sqrt(sum(M0.M.^2,1));

e3 = bsxfun(@rdivide,M0.M,S0);
% e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
e2  = [zeros(1,nMagExt); e3(3,:); -e3(2,:)];
%e2(3,~any(e2)) = 1;
e2(3,~any(abs(e2)>1e-10)) = 1;
e2  = bsxfun(@rdivide,e2,sqrt(sum(e2.^2,1)));
% e1 = e2 x e3
e1  = cross(e2,e3);

% new elliptical motion with new major and minor axes
% all phi starting phases are zero
% TODO check (+/-)
e1p = bsxfun(@times,real(Cx),e1) - bsxfun(@times,imag(Cy),e2);
e2p = bsxfun(@times,real(Cy),e2) + bsxfun(@times,imag(Cx),e1);

% check whether elliptical motion is valid
if sum((imag(Cy)./real(Cx)-imag(Cx)./real(Cy)).^2) > 1e-10
    warning('sw_spinmotion:NoElliptical','The spin precession is not elliptical!');
end

% normalize amplitudes to 1
Amp  = sqrt(sum(sum(e1p.^2+e2p.^2)))+1e-10;
e1p = e1p./Amp;
e2p = e2p./Amp;

% save results into spectra.motion field
spectra.motion.Amp1 = e1p;
spectra.motion.Amp2 = e2p;

% save additional information about the spin motion
spectra.motion.iQ = iQ;
spectra.motion.iMagnon = iMagnon;
spectra.motion.omega   = om;

end