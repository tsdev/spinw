function [r, aIdx, opInfo] = position(symOp, r0, fid, tol)
% generates symmetry equivalent positions
% 
% ### Syntax
% 
% `[r, aIdx, opInfo] = swsym.position(sym,r0)`
% 
% `[r, aIdx, opInfo] = swsym.position(sym,r0,fid)`
%
% `[r, aIdx, opInfo] = swsym.position(sym,r0,fid,tol)`
%
% ### Description
% 
% `[r, aIdx, opInfo] = swsym.position(sym, r0, fid, tol)` generates all
% symmetry equivalent atomic positions from a given space group and
% coordinates of the symmetry inequivalent atoms. If `fid` is defined, the
% result are printed onto the corresponding file.
% 
% ### Input Arguments
% 
% `sym`
% : Either the label of the space group or the index from
%   the [International Tables of Crystallography](http://it.iucr.org/A/) or
%   string containing the space group operators in the same format as used
%   in the `symmetry.dat` file (for details see [swsym.str]).
% 
% `r0`
% : Atomic position in lattice units in a matrix with dimensions of
%   $[3\times n_{atom}]$.
% 
% `fid`
% : If non-zero, the symmetry operators will be printed to the file
%   identified by `fid`, the following values are valid:
%   * `0`   no printed output (default),
%   * `1`   standard output (Command Line),
%   * `fid` text file opened before using `fid = fopen(path)`.
% 
% `tol`
% : Tolerance, distance within two atoms are considered
%   identical, default value is $10^{-5}$ lattice unit. Necessary to check
%   for badly defined atomic positions (when atoms are not exactly on the
%   symmetry element) and to avoid numerical errors.
% 
% ### Output Arguments
% 
% `r`
% : All generated atomic positions stored in a matrix with dimensions of
%   $[3\times n_{genAtom}]$.
%
% `aIdx`
% : The index of the symmetry inequivalent position for every
%   generated position, stored in a row vector with $n_{genAtom}$ number of
%   elements.
%
% `opInfo`
% : Structure with the following fields:
%   * `ismoved`     Cell, where each element contains a vector with logical
%                   values, whether the given operator moved the atom or
%                   not. Each vector has a dimensions of $[1\times n_{sym}]$, where
%                   the $n_{sym}$ is multiplicity of the general position.
%   * `opmove`      The rotation operator that moved the original atom the
%                   equivalent position stored in a matrix with dimensions
%                   of $[3\times 3\times n_{genAtom}]$.
% 
% ### See Also
% 
% [swsym.operator]
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

function r = sw_cmod(r, tol)
% modulo one with tolerance
% 
% ### Syntax
% 
% `r = sw_cmod(r, tol)`
% 
% ### Description
% 
% `r = sw_cmod(r, tol)` calculates modulo one with tolerance, numbers
% larger than $1-\epsilon$  $-\epsilon$.
% 
% ### See Also
% 
% [mod]
%

r = mod(r,1);

r(r > 1-tol) = r(r > 1-tol)-1;

end