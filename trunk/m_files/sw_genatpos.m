function [rSym, aIdx, isMoved, symName] = sw_genatpos(sym, r, fid, tol)
% [rSym, aIdx, isMoved, symName] = SW_GENATPOS(sym, r, {fid}, {tol})
% generates all symmetry equivalent atomic positions from a given symmetry
% number and coordinates of the input atoms. If print is defined, the
% result is printed onto the command window.
%
% Input:
%
% sym           Line index in the symmetry.dat file or string of the
%               symmetry operators.
% r             Atomic position in lattice units, dimensions are [3 nAtom].
% {fid}         Optional input, the file identifier to print the result.
%               To print onto the Command Window, use fid = 1; default is
%               fid = 0, no print.
% tol           Tolerance, distance within two atoms are considere, default
%               is 0.05.
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
    tol = 1e-5;
end

[symOp, symTr, symName] = sw_gencoord(sym, fid);

nAtom = size(r,2);
nSym  = size(symOp,3);
rSym  = zeros(3,0);
aIdx  = [];
isMoved  = cell(1,nAtom);

% loop over all input atoms
for ii = 1:nAtom
    % generate all equivalent atomic positions, some might overlap
    rTemp = mod(permute(sum(repmat(r(:,ii)',[3 1 nSym]).*symOp,2),[1 3 2])+symTr,1);
    % take out the overlapping positions
    isMoved{ii} = sum(bsxfun(@minus,rTemp,r(:,ii)).^2,1) > tol;
    if numel(rTemp) > 3
        [rTemp] = uniquetol(rTemp',tol,'rows')';
    end
    rSym  = [rSym rTemp]; %#ok<AGROW>
    aIdx  = [aIdx ones(1,size(rTemp,2))*ii]; %#ok<AGROW>
end

nGenAtom = numel(aIdx);

if fid
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