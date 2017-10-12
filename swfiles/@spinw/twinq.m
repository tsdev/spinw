function [Qtwin, rotQout] = twinq(obj, Q0)
% calculates equivalent Q point in twins
% 
% ### Syntax
% 
% `[qTwin, rotQ] = twinq(obj, {Q0})`
% 
% ### Description
% 
% `[qTwin, rotQ] = twinq(obj, {q0})` calculates the $Q$ values in the twin
% coordinate systems, in rlu. It also returns the rotation matrices, that
% transforms the $Q$ point from the original lattice to the selected twin
% rlu coordinate system.
% 
% ### Examples
% 
% This example Calculates the $[1,0,0]$ and $[1,1,0]$ Bragg reflections
% equivalent positions in the twins.
%
% ```
% Q1 = [1 0 0; 1 1 0];
% Q2 = cryst.twinq(Q1');
% ```
% 
% ### Input Arguments
% 
% `Q0`
% : $Q$ values in the original crystal in rlu sotred in a matrix with
% dimensions of $[3\times n_Q]$, optional.
% 
% ### Output Arguments
% 
% `Qtwin`
% : $Q$ values in the twin oordinate system in a cell element for
%           each twin.
%
% `rotQ`
% : Rotation matrices with dimensions of $[3\times 3\times n_{twin}]$.
% 
% ### See Also
% 
% [spinw] \| [spinw.addtwin]
%
% *[rlu]: Reciprocal Lattice Unit
%

if nargin == 1
    help spinw.twinq
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