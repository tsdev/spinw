function [Qtwin, rotQout] = twinq(obj, Q0)
% calculates equivalent Q point in twins
%
% [Qtwin, rotQ] = TWINQ(obj, {Q0})
%
% Qtwin are the Q values in the twin coordinate systems, in r.l.u.
% rotQ are the rotation matrices, that transforms Q points (in r.l.u.
% units) from the original crystal to the selected twin r.l.u. coordinate
% system.
%
% Input:
%
% Q0        Q values in the original crystal in r.l.u. Dimensions are
%           [3 nQ]. Optional.
%
% Output:
%
% Qtwin     Q values in the twin oordinate system in a cell element for
%           each twin.
% rotQ      Rotation matrices, dimensions are [3 3 nTwin].
%
% Example:
%
% ...
% Q1 = [1 0 0; 1 1 0];
% Q2 = cryst.twinq(Q1');
%
% This example Calculates the [1 0 0] and [1 1 0] Bragg reflections
% equivalent position in the twins.
%
% See also SPINW, SPINW.ADDTWIN.
%

if nargin == 1
    help spinw.twinq;
    return
end


% basis vectors
bv = obj.basisvector;

nTwin = size(obj.twin.vol,2);

% rotation matrices, output only if requested
rotQ = zeros(3,3,nTwin);
for ii = 1:nTwin
    %rotQ(:,:,ii) = inv(bv)*obj.twin.rotc(:,:,ii)*bv; %#ok<MINV>
    rotQ(:,:,ii) = bv\obj.twin.rotc(:,:,ii)*bv;
end
if nargout>1
    rotQout = rotQ;
end

% rotate Q points if given as input
Qtwin = cell(1,nTwin);
if nargin>1
    for ii = 1:nTwin
        Qtwin{ii} = (Q0'*rotQ(:,:,ii))';
    end
else
    Qtwin = {};
end

end