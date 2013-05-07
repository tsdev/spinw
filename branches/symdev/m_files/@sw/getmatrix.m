function aMat = getmatrix(obj, cIdx, varargin)
% gives the symmetry allowed matrices for a given coupling
%
% aMat = GETMTRIX(obj, cIdx, 'Option1', Value1, ...)
%
% Input:
%
% obj       sw object.
% cIdx      The index of the symmetry equivalent couplings.
%
% Options:
%
% fid       For printing of the allowed matrices and point group symmetries
%           onto the Command Window (fid = 1) or onto a file
%           (fid = fopen(...)). Default is fid = 0 for no printed output.
% tol       Tolerance for printing the output for the smallest matrix
%           element.
% pref      Prefactor for the symmetry allowed matrices. In the following
%           order: [A B C ... Dx Dy Dz]. First the prefactor for the
%           symmetric matrix elements, then prefactors for the
%           antisymmetric matrix elements, dimensions are [1 nSymMat].
%
% Output:
% aMat      If no prefactors are defined, aMat contains all symmetry
%           allowed elements of the coupling matrix, dimensions are 
%           [3 3 nSymMat]. If prefactor is defined, it is a single 3x3
%           matrix.
%

inpForm.fname  = {'fid' 'tol' 'pref'};
inpForm.defval = {0     1e-5  0     };
inpForm.size   = {[1 1] [1 1] [1 -1]};

param = sw_readparam(inpForm, varargin{:});
tol = param.tol;
fid = param.fid;

if nargin == 1
    help getmatrix;
    return;
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

% convert the matrices to the xyz Cartesian coordinate system
A   = obj.basisvector;
pOp = mmat(A,mmat(pOp,inv(A)));

% determine the allowed matrix elements for the first coupling
[aMat, asym] = sw_basismat(pOp,dr(:,1));

aMatS = aMat(:,:,~asym);
aMatA = aMat(:,:, asym);
dVect = permute([aMatA(2,3,:) aMatA(3,1,:) aMatA(1,2,:)],[2 3 1]);

% calculate interaction matrix using the prefactors
if any(param.pref)
    if numel(param.pref)~=size(aMat,3)
        error('sw:getmatrix:WrongInput','Wrong number of elements in pref option!');
    end
    param.pref = param.pref(:);
    aMat = sum(repmat(permute(param.pref,[2 3 1]),[3 3 1]).*aMat,3);
end

if fid
    % strings for each element in the symmetric matrix
    eStr   = cell(3);
    first  = true(3);
    firstN = true(3);
    
    for jj = 1:3
        for kk = 1:3
            eStr{jj,kk} = '';
            for ii = squeeze(find(abs(aMatS(jj,kk,:))>tol))'
                if firstN(jj,kk) && (aMatS(jj,kk,ii) > tol)
                    eStr{jj,kk} = ' ';
                    firstN(jj,kk) = false;
                elseif firstN(jj,kk)
                    firstN(jj,kk) = false;
                end
                if (aMatS(jj,kk,ii) > tol) && ~first(jj,kk)
                    eStr{jj,kk} = [eStr{jj,kk} '+']; %#ok<*AGROW>
                    first(jj,kk) = false;
                elseif first(jj,kk)
                    first(jj,kk) = false;
                end
                if abs(aMatS(jj,kk,ii)+1) < tol
                    % -1
                    eStr{jj,kk} = [eStr{jj,kk} '-' char(64+ii)];
                elseif abs(aMatS(jj,kk,ii)-1) < tol
                    % +1
                    eStr{jj,kk} = [eStr{jj,kk} char(64+ii)];
                else
                    % other number
                    eStr{jj,kk} = [eStr{jj,kk} num2str(aMatS(jj,kk,ii)) char(64+ii)];
                end
            end
            
        end
    end
    
    lStr = zeros(3);
    for ii = 1:9
        if isempty(eStr{ii}) || strcmp(eStr{ii},' ')
            eStr{ii} = ' 0';
        end
        lStr(ii) = length(eStr{ii});
    end
    % longest string
    mStr = max(lStr(:));
    % padd all the other strings to this length
    for ii = 1:9
        eStr{ii} = [eStr{ii} 32*ones(1,mStr-lStr(ii))];
    end
    
    % string of the symmetric matrix
    smatStr = '';
    for ii = 1:3
        if ii>1
            smatStr = [ smatStr '      '];
        end
        smatStr = [ smatStr '|'];
        for jj = 1:3
            smatStr = [smatStr eStr{ii,jj} '|'];
        end
        smatStr = [smatStr sprintf('\n')];
    end
    
    % strings for each element in the asymmetric matrix
    eStr   = cell(1,3);
    first  = true(1,3);
    firstN = true(1,3);
    
    for ii = 1:3
        eStr{ii} = '';
        for jj = find(abs(dVect(ii,:))>tol)
            if firstN(jj) && (dVect(ii,jj) > tol)
                eStr{ii} = ' ';
                firstN(ii) = false;
            elseif firstN(ii)
                firstN(ii) = false;
            end
            if (dVect(ii,jj) > tol) && ~first(ii)
                eStr{ii} = [eStr{ii} '+']; %#ok<*AGROW>
                first(ii) = false;
            elseif first(ii)
                first(ii) = false;
            end
            if abs(dVect(ii,jj)+1) < tol
                % -1
                eStr{ii} = [eStr{ii} '-D' char(119+jj)];
            elseif abs(dVect(ii,jj)-1) < tol
                % +1
                eStr{ii} = [eStr{ii} 'D' char(119+jj)];
            else
                % other number
                eStr{ii} = [eStr{ii} num2str(dVect(ii,jj)) 'D' char(119+jj)];
            end
        end
    end
    
    lStr = zeros(1,3);
    for ii = 1:3
        if isempty(eStr{ii}) || strcmp(eStr{ii},' ')
            eStr{ii} = ' 0';
        end
        lStr(ii) = length(eStr{ii});
    end
    % longest string
    mStr = max(lStr(:));
    % padd all the other strings to this length
    for ii = 1:3
        eStr{ii} = [eStr{ii} 32*ones(1,mStr-lStr(ii))];
    end
    
    % string of the symmetric matrix
    amatStr = ['[' eStr{1} ',' eStr{2} ',' eStr{3} ']'];
    
    % print the answer
    fprintf(fid,'\nThe symmetry analysis of the coupling between atom %d and atom %d:\n',atom1(1),atom2(1));
    fprintf(fid,' lattice translation vector: [%d,%d,%d]\n',dl(:,1));
    fprintf(fid,' distance: %5.3f Angstrom\n',norm(obj.basisvector*dr(:,1)));
    fprintf(fid,' center of bond (lattice units): [%5.3f,%5.3f,%5.3f]\n', center(:,1));
    fprintf(fid,' allowed elements in the symmetric exchange matrix:\n');
    fprintf(fid,'  S = %s\n',smatStr);
    fprintf(fid,' allowed components of the Dzyaloshinskii-Moriya vector:\n');
    fprintf(fid,'  D = %s\n\n',amatStr);
end

end