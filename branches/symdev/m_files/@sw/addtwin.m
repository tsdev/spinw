function addtwin(obj,varargin)
% adds new twins to an sw object
%
% obj = ADDTWIN(obj,'Option', Value,...)
%
% Options:
%
% axis      Defines axis of rotation to generate twins in the xyz
%           coordinate system, dimensions are [1 3].
% phi       Defines the angle of rotation to generte twins in radian,
%           several twins can be defined parallel if phi is a vector.
%           Dimensions are [1 nTwin].
% rotC      Rotation matrices, that define crystallographic twins, can be
%           given directly, dimensions are [3 3 nTwin].
% vol       Volume fractions of the twins, dimensions are [1 nTwin] default
%           value is ones(1,nTwin).
%

if nargin == 1
    help sw.addtwin;
    return;
end

inpForm.fname  = {'axis'  'phi'  'rotC'   'vol' };
inpForm.defval = {[0 0 0] 0      zeros(3) 1     };
inpForm.size   = {[1 3]   [1 -1] [3 3 -2] [1 -3]};

param = sw_readparam(inpForm, varargin{:});

% prefer axis definition over matrix
if any(param.axis)
    nTwin      = size(param.phi,2);
    param.rotC = zeros(3,3,nTwin);
    for ii = 1:nTwin
        [~, param.rotC(:,:,ii)] = sw_rot(param.axis,param.phi(ii));
    end
else
    nTwin = size(rotC,3);
end

if size(param.vol,2)<nTwin
    param.vol = ones(1,nTwin);
end

if (size(param.rotC,1)~=3) || (size(param.rotC,2)~=3)
    error('sw:addtwin:WrongInput','rotC matrix dimensions have to be [3 3 nTwin]!');
end

obj.twin.vol  = [obj.twin.vol param.vol];
obj.twin.rotc = cat(3,obj.twin.rotc, param.rotC);

end