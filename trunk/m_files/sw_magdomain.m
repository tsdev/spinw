function spectra = sw_magdomain(spectra, varargin)
% calculates the spin-spin correlation function for magnetic domains
%
% spectra = SW_MAGDOMAIN(spectra, 'Option1', Value1, 'Option2', Value2...)
%
% The function can account for magnetic domains that are related by any
% point group operation. Several domains with different volume ratios can
% be defined. The spin-spin correlation function will be rotated and summed
% according to the domains. The rotations of the magnetic domains are
% defined in the xyz coordinate system, same as the coordinate system for
% the spin-spin correlation function.
%
% Input:
%
% spectra   Calculated spin wave spectrum.
%
% Options:
%
% axis      Defines axis of rotation to generate twins in the xyz
%           coordinate system, dimensions are [1 3].
% angle     Defines the angle of rotation to generate twins in radian
%           units, several twins can be defined parallel if angle is a
%           vector. Dimensions are [1 nTwin].
% angled    Defines the angle of rotation to generate twins in degree
%           units, several twins can be defined parallel if angle is a
%           vector. Dimensions are [1 nTwin].
% rotC      Rotation matrices, that define crystallographic twins, can be
%           given directly, dimensions are [3 3 nTwin].
% vol       Volume fractions of the twins, dimensions are [1 nTwin].
%           Default value is ones(1,nTwin).
%
% Output:
%
% 'spectra' will contain the following additional/changed fields:
% Sab       The multi domain spectrum will be stored here.
% Sabraw    The original single domain spectrum is kept here, so that a
%           consecutive run of sw_magdomain will use the original single
%           domain spectrum, without the need of recalculating the full
%           spectrum.
% domVol    Volume of each domains in a vector, with dimensions of
%           [1 nDom], where nDom is the number of domains.
% domRotC   Rotation matrices for each domain, with dimensions of
%           [3 3 nDom].
%
% Example:
%
% ...
% spec = cryst.spinwave({[0 0 0] [1 0 0]});
% spec = sw_magdomain(spec,'axis',[0 0 1],'angled',[0 90 180 270]);
%
% The above example calculates the spectrum for magnetic domains that are
% related by a 90 degree rotation around the z-axis (perpendicular to the
% ab plane). All domains have equal volume.
%
% See also SW.SPINWAVE, SW.ADDTWIN, SW.TWINQ.
%

if nargin == 0
    help sw_magdomain;
    return;
end

inpForm.fname  = {'axis'  'angle' 'rotC'   'vol'  'angled'};
inpForm.defval = {[0 0 0] 0       zeros(3) 1      0       };
inpForm.size   = {[1 3]   [1 -1]  [3 3 -2] [1 -3] [1 -4]  };

param = sw_readparam(inpForm, varargin{:});

if numel(param.angled)>1 || param.angled~=0
    param.angle = param.angled *pi/180;
end

% prefer axis definition over matrix
if any(param.axis)
    nDom       = size(param.angle,2);
    param.rotC = zeros(3,3,nDom);
    for ii = 1:nDom
        [~, param.rotC(:,:,ii)] = sw_rot(param.axis,param.angle(ii));
    end
else
    nDom = size(param.rotC,3);
end

if size(param.vol,2)<nDom
    param.vol = ones(1,nDom);
end

if (size(param.rotC,1)~=3) || (size(param.rotC,2)~=3)
    error('sw_magdomain:WrongInput','rotC matrix dimensions have to be [3 3 nDom]!');
end

% param.vol, param.rotC defined
% add together the domains
vol = param.vol/sum(param.vol);
rotC = param.rotC;

SabDom = spectra.Sab*0;

if isfield(spectra,'Sabraw')
    for ii = 1:nDom
        SabDom = SabDom + mmat(mmat(rotC(:,:,ii),spectra.Sabraw),rotC(:,:,ii)')*vol(ii);
    end
else
    for ii = 1:nDom
        SabDom = SabDom + mmat(mmat(rotC(:,:,ii),spectra.Sab),rotC(:,:,ii)')*vol(ii);
    end
end

% save the new matrix
spectra.Sabraw = spectra.Sab;
spectra.Sab    = SabDom;

% save the magnetic domain information
spectra.domVol = vol;
spectra.domRotC = rotC;

end