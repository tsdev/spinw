function spectra = sw_conv(spectra, varargin)
% spectra = SW_CONV(spectra, 'Option1', Value1, ...) convolutes different
% magnetic scattering cross sections with the spin wave dispersion.
%
% Input:
%
% spectra   Input structure, contains calculated correlation functions and
%           different cross sections.
%
% Options:
%
% convmode  Selects the previously calculated intensity component to be
%           convoluted. The possible options are:
%               'Sperp' convolutes the magnetic neutron scattering
%                       intensity (<Sperp * Sperp> expectation value).
%                       Default.
%               'Sab'   convolutes the selected components of the spin-spin
%                       correlation function. Letter a and b can be 'x',
%                       'y' or 'z'. For example: 'Sxx' will convolute the
%                       xx component of the correlation function with the
%                       dispersion. xyz is the standard coordinate system,
%                       see online documentation of sw.
%               'Mab'   convolutes the selected components of the spin-spin
%                       correlation function. Letter a and b can be 'x',
%                       'y' or 'z'. For example: 'Sxx' will convolute the
%                       xx component of the correlation function with the
%                       dispersion. The xyz coordinates are in the
%                       Blume-Maleev coordinate system, see below.
%               'Pab'   convolutes the selected element of the polarisation
%                       matrix. Letter a and b can be 'x', 'y' or 'z'. For
%                       example: 'Pyy' will convolute the yy component of
%                       the polarisation matrix with the dispersion. The
%                       xyz coordinates are in the Blume-Maleev coordinate
%                       system, see below.
%               'Pa'    convolutes the intensity of the simulated polarised
%                       neutron scattering, with inciden polarisation of
%                       Pa. Letter a can be 'x', 'y' or 'z'. For example:
%                       'Py' will convolute the scattering intensity
%                       simulated for incident polarisation Pi || y. The
%                       xyz coordinates are in the Blume-Maleev coordinate
%                       system, see below.
%               'fName' where fName is one of the field names of the input
%                       structure spectra. This field should contain a
%                       matrix with size [nMode nHkl].
%           Any linear combination of the above are allowed, for example:
%           'Sxx+2*Syy' convolutes the linear combination of the xx
%           component of the spin-spin correlation function and the yy
%           component.
%           Several cross section can be convoluted and stored
%           independently, if convmode is a cell array containing strings
%           each containing any linear combination of cross sections as
%           above, the cell array needs to have size [1 nCell].
%
% Evect     Vector, defined the center of the energy bins of the convoluted
%           output, size IS [1 nE]. The energy units are defined by the
%           unit.kB property of the sw object. Default is
%           linspace(0,1,100).
% T         Temperature to calculate the Bose factor in units
%           depending on the Boltzmann constant (sw.unit.kB). Default is
%           zero (no Bose factor).
% sumtwin   If true, the spectra of the different twins will be summed
%           together. Default is true.
%
% The Blume-Maleev coordinate system is a cartesian coordinate system
% with (xBM, yBM and zBM) basis vectors as follows:
%           xBM    parallel to the momentum transfer Q,
%           yBM    perpendicular to xBM in the scattering plane,
%           zBM    perpendicular to the scattering plane.
%
%
% Output:
%
% spectra contains the following additional fields beside the input:
%
% swConv    Stores the selected cross section convoluted with the
%           dispersion, size is [nE nHkl].
% swInt     Stores the selected cross section, size is [nMode nHkl].
%
% If convmode is a cell array or the spectra of multiple domains are
% convoluted separately, swConv and swInt will be packaged into a cell. The
% dimensions of the cell are [nConv nTwin].
%
% T         Input temperature.
% convmode  Input convolution mode.
% Evect     Input energy bin.
%
% See also SW.SPINWAVE, SW.SWINC, SW_NEUTRON.
%

if nargin == 0
    help sw_conv;
    return;
end

inpForm.fname  = {'Evect'           'T'   'convmode' 'sumtwin'};
inpForm.defval = {linspace(0,1,100) 0     'Sperp'    true     };
inpForm.size   = {[1 -1]            [1 1] [1 -2]     [1 1]    };

param = sw_readparam(inpForm, varargin{:});

% parse the convmode string
if iscell(param.convmode)
    nConv = numel(param.convmode);
    parsed = cell(1,nConv);
    for ii = 1:numel(param.convmode)
        parsed{ii} = sw_parstr(param.convmode{ii});
    end
elseif isstruct(param.convmode)
    nConv  = 1;
    parsed = {param.convmode};
    param.convmode = {parsed{1}.string};
else
    nConv = 1;
    parsed = {sw_parstr(param.convmode)};
    param.convmode = {param.convmode};
end

% pack all cross section into a cell for easier looping
if iscell(spectra.omega)
    nTwin = numel(spectra.omega);
    omega = spectra.omega;
    Sab   = spectra.Sab;
    intP  = spectra.intP;
    Pab   = spectra.Pab;
    Mab   = spectra.Mab;
    Sperp = spectra.Sperp;

else
    nTwin = 1;
    omega = {spectra.omega};
    Sab   = {spectra.Sab};
    intP  = {spectra.intP};
    Pab   = {spectra.Pab};
    Mab   = {spectra.Mab};
    Sperp = {spectra.Sperp};
    
end

% extract the requested cross section
nMode = size(omega{1},1);
nHkl  = size(omega{1},2);
% DSF stores the intensity that is going to be convoluted
DSF = cell(nConv,nTwin);

for tt = 1:nTwin
    for ii = 1:numel(parsed)
        par0 = parsed{ii};
        DSF{ii,tt} = zeros(nMode, nHkl);
        for jj = 1:length(par0.type)
            switch par0.type{jj}(1)
                case 1
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*Sperp{tt};
                case 2
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*permute(Sab{tt}(par0.type{jj}(2),par0.type{jj}(3),:,:),[3 4 1 2]);
                case 3
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*permute(Mab{tt}(par0.type{jj}(2),par0.type{jj}(3),:,:),[3 4 1 2]);
                case 4
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*permute(Pab{tt}(par0.type{jj}(2),par0.type{jj}(3),:,:),[3 4 1 2]);
                case 5
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*permute(intP{tt}(par0.type{jj}(2),:,:),[2 3 1]);
                otherwise
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*spectra.(par0.type{jj});
            end
        end
    end
end

% Create vector for energy values, and put extra value below minimum and
% above maximum for easy indexing swConv.
Evect     = sort(param.Evect);
epsilon   = 1e-8;

% Separate the negative and positive energies.
EvectNeg  = Evect(Evect<0);
if ~isempty(EvectNeg)
    EvectNeg = [EvectNeg(1)-epsilon; EvectNeg(:); EvectNeg(end)+epsilon];
else
    EvectNeg = [-epsilon; epsilon];
end
nENeg = length(EvectNeg);

EvectPos  = Evect(Evect>=0);
if ~isempty(EvectPos)
    EvectPos = [EvectPos(1)-epsilon; EvectPos(:); EvectPos(end)+epsilon];
else
    EvectPos = [-epsilon; epsilon];
end
nEPos = length(EvectPos);

% Create indices in the matrix by searching for the closest value, size
% nMode x nHkl. Put all the modes to the positive side for magnon creation.
% The negative side will be the same, however with different Bose factor
% for non-zero temperature.
idxNeg = cell(1,nTwin);
idxPos = cell(1,nTwin);
for tt = 1:nTwin
    [~, idxNeg{tt}] = min(abs(repmat(-abs(omega{tt}),[1 1 nENeg])-repmat(permute(EvectNeg,[2 3 1]),[nMode nHkl 1])),[],3);
    [~, idxPos{tt}] = min(abs(repmat( abs(omega{tt}),[1 1 nEPos])-repmat(permute(EvectPos,[2 3 1]),[nMode nHkl 1])),[],3);
    
    % Creates indices in the swConv matrix.
    idxNeg{tt} = idxNeg{tt} + repmat((0:nHkl-1).*nENeg,[nMode 1]);
    idxNeg{tt} = idxNeg{tt}(:);
    idxPos{tt} = idxPos{tt} + repmat((0:nHkl-1).*nEPos,[nMode 1]);
    idxPos{tt} = idxPos{tt}(:);
end

% Sums up the intensities in DSF into swConv.
swConvNeg = cell(nConv,nTwin);
swConvPos = cell(nConv,nTwin);

for tt = 1:nTwin
    for ii = 1:nConv
        swConvNeg{ii,tt} = reshape(accumarray(idxNeg{tt},DSF{ii,tt}(:),[nENeg*nHkl 1]),[nENeg nHkl]);
        swConvPos{ii,tt} = reshape(accumarray(idxPos{tt},DSF{ii,tt}(:),[nEPos*nHkl 1]),[nEPos nHkl]);
    end
end

% Bose temperature factor for magnons
if param.T==0
    nNeg = 1+0*EvectNeg(2:(end-1));
else
    nNeg = 1./(exp(abs(EvectNeg(2:(end-1)))/(spectra.obj.unit.kB*param.T))-1);
end

if param.T == 0
    nPos = 1+0*EvectPos(2:(end-1));
else
    nPos = 1./(exp((abs(EvectPos(2:(end-1)))+1e-8)/(obj.unit.kB*param.T))-1)+1;
end

% Unite the negative and positive energy transfer sides.
swConv = cell(nConv,nTwin);
for tt = 1:nTwin
    for ii = 1:nConv
        swConv{ii,tt} = [bsxfun(@times,swConvNeg{ii,tt}(2:(end-1),:),nNeg); bsxfun(@times,swConvPos{ii,tt}(2:(end-1),:),nPos)];
    end
end

% sum up twin spectra if requested
if param.sumtwin
    for tt = 1:nTwin
        for ii = 1:nConv
            swConv{ii,1} = swConv{ii,1} + swConv{ii,tt};
            DSF{ii,1}    = DSF{ii,1} + DSF{ii,tt};
        end
    end
    swConv = swConv(:,1);
    DSF    = DSF(:,1);
end


spectra.T        = param.T;
spectra.Evect    = param.Evect;

if numel(swConv) == 1
    spectra.swConv   = swConv{1};
    spectra.swInt    = DSF{1};
    spectra.convmode = param.convmode{1};
else
    spectra.swConv   = swConv;
    spectra.swInt    = DSF;
    spectra.convmode = param.convmode;
end

end