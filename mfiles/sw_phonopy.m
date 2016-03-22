function spectra = sw_phonopy(fName)
% extract dispersion from Phonopy output yaml file
%
% spectra = SW_PHONOPY(fName)
%
% Energies are converted assuming Phonopy output is THz to meV. Since
% Phonopy uses different Descartes coordinate system, an automatic
% conversion to SpinW coordinate system is made. Also a fictious magnetic
% structure is setup to set the nMagExt value equal to the number of phonon
% modes divided by two. As a side effect, the "q-position" field in the
% .yaml file is renamed to "qposition".
%
% Input:
%
% fName     String that contains the location of the Phonopy output .yaml
%           file.
%
% Output:
%
% spectra   Spectrum that is compatible with all the standard SpinW
%           functions.
%
% See also SW_XRAY, SW_NEUTRON, SW_EGRID, SW_INSTRUMENT, SW_PLOTSPEC.
%

if nargin == 0
    help sw_phonopy
    return
end

% change all q-points --> qpoints before importing
system(['perl -i -pe ''s/q-position/qposition/g'' ' fName]);

% load the file using Yaml
fprintf('Loading the .yaml file...\n')
band = YAML.read(fName);
fprintf('Finished.\n')

nHkl  = numel(band.phonon);
nMode = numel(band.phonon(1).band);

spectra.hkl   = zeros(3,nHkl);
spectra.omega = zeros(nMode,nHkl);
spectra.Sab   = zeros(band.natom*3,nMode,nHkl);

for ii = 1:nHkl
    spectra.hkl(:,ii) = band.phonon(ii).qposition';
    for jj = 1:nMode
        spectra.omega(jj,ii) = band.phonon(ii).band(jj).frequency;
        spectra.Sab(:,jj,ii) = sum(bsxfun(@times,cell2mat(band.phonon(ii).band(jj).eigenvector),[1 1i]),2);
    end
end

% reshape the modes
spectra.Sab = reshape(spectra.Sab,3,[],nMode,nHkl);

spectra.omega = spectra.omega*sw_converter(1,'THz','meV');

spectra.obj = spinw;
% generate the object from the basis vectors
R = spectra.obj.genlattice('bv',inv(band.reciprocal_lattice));
% add atoms
for ii = 1:band.natom
    spectra.obj.addatom('r',band.atomReducedPositions(ii,:),'label',band.atomNames{ii},'S',1);
end

% fake magnetic structure to generate the right nMagExt = nAtom value
spectra.obj.mag_str.S = zeros(3,nMode/2);
% rotate the components of the polarisation using R
spectra.Sab = permute(mmat(R,permute(spectra.Sab,[1 5 2:4])),[1 3:5 2]);

spectra.norm = false;
% Angstrom^-1 units for Q
spectra.hklA = 2*pi*(spectra.hkl'/spectra.obj.basisvector)';

end