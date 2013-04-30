function [rSym, symName] = sw_genatpos(sym, r, fid, tol)
% [rSym, symName] = SW_GENATPOS(sym, r, {fid}, {tol}) generates all symmetry
% equivalent atomic positions from a given symmetry number and coordinates
% of the input atoms. If print is defined, the result is printed onto
% the command window.
%
% Input:
%
% sym           Line index in the symmetry.dat file or string of the
%               symmetry operators.
% r             Atomic position in lattice units, dimensions are [3 nAtom].
% {fid}         Optional input, the file identifier to print the result.
%               To print onto the Command Window, use fid = 1; default is
%               fid = 0, no print.
%
% Output:
%
% rSym          Symmetry equivalent atomic positions, dimensions are
%               [3 nGenAtom] if the input was a single atom. If more atoms
%               are defined, the output is packaged into a cell, with
%               dimensions of [1 nAtom]. Every element of the cell contains
%               a matrix with the symmetry equivalent positions.
% symName       String, the name of  the space group.
% tol           Tolerance, distance within two atoms are considere, default is 0.05.
%
% See also SW, SW.ATOM, SW.MATOM, SW_GENCOUPLING, SW_GENCOORD, SW_GENSYM.
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
rSym  = cell(1,0);

% loop over all input atoms
for ii = 1:nAtom
    % generate all equivalent atomic positions, some might overlap
    rTemp = mod(permute(sum(repmat(r(:,ii)',[3 1 nSym]).*symOp,2),[1 3 2])+symTr,1);
    % take out the overlapping positions
    rSym{ii} = consolidator(rTemp',[],[],tol)';
end

if fid
    for ii = 1:nAtom
        fprintf(fid,'\nAtomic coordinates generated for: (%5.3f %5.3f %5.3f)\n',r(:,ii));
        for jj = 1:size(rSym{ii},2)
            fprintf(fid,'r%i (%5.3f %5.3f %5.3f)\n',jj,rSym{ii}(:,jj));
        end
    end
    fprintf(fid,'\n');
end

if nAtom == 1
    rSym = rSym{1};
end

end