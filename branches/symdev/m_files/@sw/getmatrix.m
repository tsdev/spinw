function [aMat] = getmatrix(obj, cIdx, fid)
% gives the symmetry allowed matrices for a given coupling
%
% [aMat] = GETMTRIX(obj, cIdx, {fid})
%
% Input:
%
% obj       sw object.
% cIdx      The index of the symmetry equivalent couplings.
% {fid}     For printing of the allowed matrices and point group symmetries
%           onto the Command Window (fid = 1) or onto a file
%           (fid = fopen(...)). Default is fid = 0 for no printed output.
%
tol = 1e-5;

if nargin == 1
    help getmatrix;
    return;
end

if nargin < 3
    fid = 0;
end

iSel = obj.coupling.idx == cIdx;

mAtom = obj.matom;
% indices of atoms in selected couplings
atom1  = obj.coupling.atom1(iSel);
atom2  = obj.coupling.atom2(iSel);
% positions of atoms in selected couplings
r1     = mAtom.r(:,atom1);
r2     = mAtom.r(:,atom2);
% lattice translation vector between interacting atoms
dl     = double(obj.coupling.dl(:,iSel));
% vector pointing from atom1 to atom2
dr     = r2 + dl - r1;
% centers of the selected couplings
center = mod((r1+r2+dl)/2,1);
% get the point group symmetry operators of the center of the first
% coupling
pOp = sw_pointsym(obj.lattice.sym,center(:,1));

% determine the allowed matrix elements for the first coupling
aMat = sw_basismat(pOp,dr(:,1));

if fid
    % print output
    naMat = size(aMat,3);
    
    % strings for each element in the matrix
    eStr = cell(3);
    lStr = zeros(3);
    for jj = 1:3
        for kk = 1:3
            first = true;
            firstN = true;
            for ii = 1:naMat
                if abs(aMat(jj,kk,ii)) > tol
                    if firstN && (aMat(jj,kk,ii) > 0)
                        eStr{jj,kk} = ' ';
                        firstN = false;
                    elseif firstN
                        firstN = false;
                    else
                        eStr{jj,kk} = '';
                    end
                    if (aMat(jj,kk,ii) > 0) && ~first
                        eStr{jj,kk} = [eStr{jj,kk} '+']; %#ok<*AGROW>
                        first = false;
                    end
                    if abs(aMat(jj,kk,ii)+1) < tol
                        % -1
                        eStr{jj,kk} = [eStr{jj,kk} '-' char(64+ii)];
                    elseif abs(aMat(jj,kk,ii)-1) < tol
                        % +1
                        eStr{jj,kk} = [eStr{jj,kk} char(64+ii)];
                    else
                        % other number
                        eStr{jj,kk} = [eStr{jj,kk} num2str(aMat(jj,kk,ii)) char(64+ii)];
                    end
                end
            end
            if isempty(eStr{jj,kk})
                eStr{jj,kk} = ' 0';
            end
            lStr(jj,kk) = length(eStr{jj,kk});
        end
    end
    % longest string
    mStr = max(lStr(:));
    % padd all the other strings to this length
    for ii = 1:9
        eStr{ii} = [eStr{ii} 32*ones(1,mStr-lStr(ii))];
    end
    
    % print the matrix
    for ii = 1:3
        fprintf(fid,'|');
        for jj = 1:3
            fprintf(fid,'%s|',eStr{ii,jj});
        end
        fprintf(fid,'\n');
    end
end


end










