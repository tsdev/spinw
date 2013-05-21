function [rSym, aIdx, isMoved, opMove, symName] = sw_genatpos(sym, r, fid, tol)
% [rSym, aIdx, isMoved, opMove, symName] = SW_GENATPOS(sym, r, fid, tol)
% generates all symmetry equivalent atomic positions from a given symmetry
% number and coordinates of the input atoms. If print is defined, the
% result is printed onto the command window.
%
% Input:
%
% sym           Line index in the symmetry.dat file or string of the
%               symmetry operators or a cell containing the output of
%               sw_gencoord().
% r             Atomic position in lattice units, dimensions are [3 nAtom].
% fid           Optional input, the file identifier to print the result.
%                   0   No output printed (Default)
%                   1   Output printed onto the screen (Command Window)
%                   fid Use with the following command: fid = fopen(path)
% tol           Tolerance, distance within two atoms are considered
%               identical, default is 0.05 Angstrom. Necessary for badly
%               defines atomic positions (when atoms are not exactly on the
%               symmetry element) and to avoid numerical errors.
%
% Output:
%
% rSym          All generated atomic positions, dimensions are
%               [3 nGenAtom].
% aIdx          The index of the symmetry inequivalent position for every
%               generated position, dimensions are [1 nGenAtom].
% isMoved       Cell each contains a vector with logical value, whether the
%               given operator moved the atom or not. Each vector has a
%               dimensions of [1 nSym], where the nSym is multiplicity of
%               the general position.
% opMove        The rotation operators for every moved atom, dimensions are
%               [3 3 nGenAtom].
% symName       String, the name of  the space group.
%
% See also SW, SW.ATOM, SW.MATOM, SW.GENCOUPLING, SW_GENCOORD, SW_GENSYM, SW_POINTSYM.
%

if nargin == 0
    help sw_genatpos;
    return;
end

if nargin == 2
    fid = 0;
end

if nargin < 4
    tol = 0.05;
end

if ~iscell(sym)
    [symOp, symTr, symName] = sw_gencoord(sym, fid);
else
    symOp = sym{1};
    symTr = sym{2};
end

nAtom  = size(r,2);
nSym   = size(symOp,3);
rSym   = zeros(3,0);
aIdx   = [];
isMoved  = cell(1,nAtom);
opMove = zeros(3,3,0);

% loop over all input atoms
for ii = 1:nAtom
    % generate all equivalent atomic positions, some might overlap
    rTemp  = mod(permute(sum(repmat(r(:,ii)',[3 1 nSym]).*symOp,2),[1 3 2])+symTr,1);
    % take out the overlapping positions, using modulo with tolerance
    isMoved{ii} = sum(bsxfun(@minus,sw_cmod(rTemp,tol),sw_cmod(r(:,ii),tol)).^2,1) > tol^2;
    if numel(rTemp) > 3
        % select unique atomic positions and the indices
        [rTemp, idxF] = sw_uniquetol(sw_cmod(rTemp,tol),tol);
    else
        idxF = 1;
    end
    rSym   = [rSym rTemp]; %#ok<AGROW>
    aIdx   = [aIdx ones(1,size(rTemp,2))*ii]; %#ok<AGROW>
    opMove = cat(3, opMove, symOp(:,:,idxF));
end

nGenAtom = numel(aIdx);

if fid ~= 0
    idx = 1;
    for ii = 1:nAtom
        fprintf(fid,'\nAtomic coordinates generated for: (%5.3f %5.3f %5.3f)\n',r(:,ii));
        idx2 = 1;
        while (aIdx(idx) == ii) && (idx<nGenAtom)
            fprintf(fid,'r%i (%5.3f %5.3f %5.3f)\n',idx2,rSym(:,idx));
            idx2 = idx2 + 1;
            idx  = idx  + 1;
        end
    end
    fprintf(fid,'\n');
end

end