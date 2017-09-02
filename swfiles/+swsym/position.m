function [r, aIdx, opInfo] = position(symOp, r0, fid, tol)
% generates symmetry equivalent positions
%
% [r, aIdx, opInfo] = SWSYM.POSITION(sym, r0, fid, tol)
% 
% It generates all symmetry equivalent atomic positions from a given
% symmetry number and coordinates of the input atoms. If fid is defined,
% the result is printed onto the command window.
%
% Input:
%
% symOp         Matrix containing the symmetry operators 
% r0            Atomic position in lattice units, dimensions are [3 nAtom].
% fid           Optional input, the file identifier to print the result.
%                   0   No output printed (Default)
%                   1   Output printed onto the screen (Command Window)
%                   fid Use with the following command: fid = fopen(path)
% tol           Tolerance, distance within two atoms are considered
%               identical, default is 1e-5 lattice unit. Necessary for
%               badly defined atomic positions (when atoms are not exactly
%               on the symmetry element) and to avoid numerical errors.
%
% Output:
%
% rSym          All generated atomic positions, dimensions are
%               [3 nGenAtom].
% aIdx          The index of the symmetry inequivalent position for every
%               generated position, dimensions are [1 nGenAtom].
% opInfo        Structure with the following fields:
%   ismoved         Cell, where each element contains a vector with logical
%                   value, whether the given operator moved the atom or
%                   not. Each vector has a dimensions of [1 nSym], where
%                   the nSym is multiplicity of the general position.
%   opmove          The rotation operators for every moved atom, dimensions
%                   are [3 3 nGenAtom].
%
% See also SPINW, SWSYM.OPERATOR, SPINW.ATOM, SPINW.MATOM.
%

if nargin == 0
    help swsym.position
    return
end

% for empty operator, just return the original positions
if isempty(symOp)
    symOp = [eye(3) zeros(3,1)];
end

if nargin == 2
    fid = 0;
end

if nargin < 4
    tol = 1e-5;
end

if size(r0,1)~=3
    error('position:WrongInput','The positions have to be in column vector format!')
end

nAtom   = size(r0,2);
r    = zeros(3,0);
aIdx    = [];
isMoved = cell(1,nAtom);
opMove  = zeros(3,3,0);

% loop over all input atoms
for ii = 1:nAtom
    % generate all equivalent atomic positions, some might overlap
    % corrected for R*pos+T
    rTemp  = permute(mod(mmat(symOp(:,1:3,:),r0(:,ii))+symOp(:,4,:),1),[1 3 2]);
    % take out the overlapping positions, using modulo with tolerance
    isMoved{ii} = sum(bsxfun(@minus,sw_cmod(rTemp,tol),sw_cmod(r0(:,ii),tol)).^2,1) > tol^2;
    if numel(rTemp) > 3
        % select unique atomic positions and the indices
        [rTemp, idxF] = sw_uniquetol(sw_cmod(rTemp,tol),tol);
    else
        idxF = 1;
    end
    r   = [r rTemp]; %#ok<AGROW>
    aIdx   = [aIdx ones(1,size(rTemp,2))*ii]; %#ok<AGROW>
    opMove = cat(3, opMove, symOp(:,1:3,idxF));
end

opInfo.ismoved = isMoved;
opInfo.opmove  = opMove;

nGenAtom = numel(aIdx);

if fid ~= 0
    idx = 1;
    for ii = 1:nAtom
        fprintf(fid,'\nAtomic coordinates generated for: (%5.3f %5.3f %5.3f)\n',r0(:,ii));
        idx2 = 1;
        while idx<=nGenAtom && aIdx(idx) == ii
            fprintf(fid,'R%02i = (%5.3f %5.3f %5.3f)\n',idx2,r(:,idx));
            idx2 = idx2 + 1;
            idx  = idx  + 1;
        end
    end
end

end