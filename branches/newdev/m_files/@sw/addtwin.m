function addtwin(obj,varargin)
% adds new twins to an sw object
%
% ADDTWIN(obj,'Option', Value,...)
%
% Input:
%
% obj       sw class object.
%
% Options:
%
% axis      Defines axis of rotation to generate twins in the xyz
%           coordinate system, dimensions are [1 3].
% phi       Defines the angle of rotation to generate twins in radian
%           units, several twins can be defined parallel if phi is a
%           vector. Dimensions are [1 nTwin].
% phid      Defines the angle of rotation to generate twins in degree
%           units, several twins can be defined parallel if phi is a
%           vector. Dimensions are [1 nTwin].
% rotC      Rotation matrices, that define crystallographic twins, can be
%           given directly, dimensions are [3 3 nTwin].
% vol       Volume fractions of the twins, dimensions are [1 nTwin].
%           Default value is ones(1,nTwin).
%
% Output:
%
% The function adds extra entries in the 'twin' field of the obj sw object.
%
% Example:
%
% ...
% cryst.addtwin('axis',[0 0 1],'phid',[60 120],'vol',[1 1])
%
% This will add two extra crystallographic twins to the crystal, with the
% original orientation there are will be three twins with equal volumes.
%
% See also SW, SW.TWINQ.
%

if nargin == 1
    help sw.addtwin;
    return;
end

inpForm.fname  = {'axis'  'phi'  'rotC'   'vol'  'phid'};
inpForm.defval = {[0 0 0] 0      zeros(3) 1      0     };
inpForm.size   = {[1 3]   [1 -1] [3 3 -2] [1 -3] [1 -4]};

param = sw_readparam(inpForm, varargin{:});

if numel(param.phid)>1 || param.phid~=0
    param.phi = param.phid *pi/180;
end

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