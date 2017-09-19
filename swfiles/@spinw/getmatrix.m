function [aMatOut, paramOut, pOpOut] = getmatrix(obj, varargin)
% determines the symmetry allowed tensor elements
% 
% ### Syntax
% 
% `amat = getmatrix(obj,Name,Value)`
% 
% ### Description
% 
% `amat = getmatrix(obj,Name,Value)` determines the symmetry allowed
% elements of the exchange, single-ion anistropy and g-tensor. For bonds,
% the code first determines the point group symmetry on the center of the
% bond and calculates the allowed eelements of the exchange tensor
% accordingly. For anisotropy and g-tensor, the point group symmetry of the
% selected atom is considered. For example the code can correctly calculate
% the allowed Dzyaloshinskii-Moriya vectors.
% 
% ### Examples
% 
% To following code will determine the allowed anisotropy matrix elements
% in the $C4$ point group (the symmetry at the $(0,0,0)$ atomic position).
% The allowed matrix elements will be `diag([A A B])`:
%
% ```
% >>cryst = spinw;
% >>cryst.genlattice('sym','P 4')
% >>cryst.addatom('r',[0 0 0],'label','MCu2')
% >>cryst.addmatrix('label','A','value',1)
% >>cryst.gencoupling
% >>cryst.addaniso('A')
% >>cryst.getmatrix('mat','A');
% ```
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Name-Value Pair Arguments
%
% At least one of the following option has to be defined:
% 
% `mat`
% : Label or index of a matrix that is already assigned to
%   a bond, anisotropy or g-tensor, e.g. `J1`.
% 
% `bond`
% : Index of the bond in `spinw.coupling.idx`, e.g. 1 for first neighbor
%   bonds.
% 
% `aniso`
% : Label or index of the magnetic atom that has a single ion
%   anisotropy matrix is assigned, e.g. `Cr1` if `Cr1` is a magnetic atom.
% 
% `gtensor`
% : Label or index of the magnetic atom that has a g-tensor is 
%   assigned.
%
% Optional inputs:
% 
% `subIdx`
% : Selects a certain bond, within equivalent bonds. Default value is 1.
%
% `tol`
% : Tolerance for printing the output for the smallest matrix
%   element.
% 
% `pref`
% : If defined `amat` will contain a single $[3\times 3]$ matrix by
%   multuplying the calculated tensor components with the given prefactors.
%   Thus `pref` should contain the same number of elements as the number of
%   symmetry allowed tensor components. Alternatively, if only a few of the
%   symmetry allowed matrices have non-zero prefactors, use e.g. 
%   `{[6 0.1 5 0.25]}` which means, the 6th symmetry allowed matrix have
%   prefactor 0.1, the 5th symmetry allowed matrix have prefactor 0.25.
%   Since Heisenberg isotropic couplings are always allowed, a cell with a
%   single element will create a Heisenberg coupling, e.g. `{0.1}, which is
%   identical to `obj.matrix.mat = eye(3)*0.1`. For Dzyaloshinskii-Moriya
%   interactions (antisymmetric exchange matrices), use a three element
%   vector in a cell, e.g. `pref = {[D1 D2 D3]}`. In this case, these will
%   be the prefactors of the 3 antisymmetric allowed matrices. In
%   case no crystal symmetry is defined, these will define directly the
%   components of the  Dzyaloshinskii-Moriya interaction in the xyz
%   coordinate system.
%
%   {{note Be carefull with the sign of the Dzyaloshinskii-Moriya
%   interaction, it depends on the counting order of the two interacting
%   atoms! For single-ion anisotropy and g-tensor antisymmetric matrices
%   are forbidden in any symmetry.}}
% 
% ### Output Arguments
% 
% `aMat`
% : If no prefactors are defined, `aMat` contains all symmetry
%   allowed elements of the selected tensor, dimensions are $[3\times 3\times n_{symmat}]$.
%   If a prefactor is defined, it is a single $[3\times 3]$ matrix, that is
%   a sum of all symmetry allowed elemenets multiplied by the given
%   prefactors.
% 
% ### See Also
% 
% [spinw.setmatrix]
%

inpForm.fname  = {'mat'      'aniso' 'bond' 'tol' 'pref' 'gtensor' 'subIdx' };
inpForm.defval = {zeros(1,0) 0       0       1e-5  []     0        1        };
inpForm.size   = {[1 -1]     [1 -2]  [1 1]   [1 1] [1 -2] [1 -3]   [1 1]    };
inpForm.soft   = {false      false   false   false true   false    false    };

param0 = sw_readparam(inpForm, varargin{:});
param  = param0;

tol = param.tol;

if nargin == 1
    help spinw.getmatrix
    return
end

if ~obj.symmetry
    % there are no symmetry operators to analyse
    if nargout>0
        aMatOut  = zeros(3,3,0);
        paramOut = [];
        pOpOut   = zeros(3,3,0);
    end
    warning('spinw:getmatrix:NoSymmetry','The SpinW object does not contain symmetry!')
    return
end
    
bondIdx = param.bond;

if ischar(param.aniso)
    % find atom labels
    anisoIdx = find(cellfun(@(C)~isempty(C),strfind(obj.unit_cell.label, param.aniso)));
    if isempty(anisoIdx)
        anisoIdx = 0;
    elseif numel(anisoIdx) == 1
        anisoIdx = anisoIdx(1);
    else
        error('spinw:getmatrix:WrongInput','The given atom label is not unique!')
    end
else
    anisoIdx = param.aniso;
end
    
if ischar(param.gtensor)
    % find atom labels
    gIdx = find(cellfun(@(C)~isempty(C),strfind(obj.unit_cell.label, param.gtensor)));
    if isempty(gIdx)
        gIdx = 0;
    elseif numel(gIdx)==1
        gIdx = gIdx(1);
    else
        error('spinw:getmatrix:WrongInput','The given atom label is not unique!')
    end
else
    gIdx = param.gtensor;
end
    
if ischar(param.mat) 
    matIdx = find(strcmp(obj.matrix.label, param.mat));
    if isempty(matIdx)
        error('spinw:getmatrix:WrongInput','Matrix label cannot be found (case sensitive)!');
    elseif numel(matIdx) > 1
        error('spinw:getmatrix:WrongInput','Multiple identical matrix labels exist!');
    end
else
    matIdx = param.mat;
end

% Identify the coupling, anisotropy or g-tensor
if matIdx ~= 0
    
    % search to which coupling/anisotropy the matrix is assigned to
    [~, selCpIdx] = find(obj.coupling.mat_idx == matIdx);
    [~, selAnIdx] = find(obj.single_ion.aniso == matIdx);
    [~, selGIdx ] = find(obj.single_ion.g     == matIdx);
    
    Cpidx = unique(obj.coupling.idx(selCpIdx));
    Anidx = unique(obj.matom.idx(selAnIdx));
    Gidx  = unique(obj.matom.idx(selGIdx));
    
    sumNum = numel(Cpidx) + numel(Anidx) + numel(Gidx);
    
    if isempty(Cpidx) && isempty(Anidx) && isempty(Gidx)
        error('spinw:getmatrix:WrongInput','Matrix is not assigned to any coupling/anisotropy/g-tensor!');
    elseif sumNum > 1
        error('spinw:getmatrix:WrongInput','Matrix is assigned to multiple coupling/anisotropy/g-tensor!');
    else
        if ~isempty(Cpidx)
            bondIdx = Cpidx;
        elseif ~isempty(Anidx)
            anisoIdx = Anidx;
        else
            gIdx = Gidx;
        end
    end
end

if bondIdx ~= 0
    % Coupling is defined
    iSel = find(obj.coupling.idx == bondIdx);
    
    if isempty(iSel)
        error('spinw:getmatrix:WrongInput','The given bond index does not existst, check your input!');
    end

    if param.subIdx ~=1
        iSel = iSel(param.subIdx);
    end

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
    matIdx0 = obj.coupling.mat_idx(:,iSel(1));
    if sum(matIdx0>0) == 1
        matIdx = sum(matIdx0);
    end
elseif anisoIdx ~= 0
    % Anisotropy is defined
    center = obj.unit_cell.r(:,anisoIdx);
    dr     = 0;
    
    matIdx0 = obj.single_ion.aniso(obj.matom.idx == anisoIdx);
    if isempty(matIdx0)
        error('spinw:getmatrix:WrongInput','The given obj.matom.idx does not exists, check your input!');
    end
    if matIdx0(1) > 0
        matIdx = matIdx0(1);
    end
else
    % g-tensor is defined
    center = obj.unit_cell.r(:,gIdx);
    dr     = 0;
    
    matIdx0 = obj.single_ion.g(obj.matom.idx == gIdx);
    if isempty(matIdx0)
        error('spinw:getmatrix:WrongInput','The given obj.matom.idx does not exists, check your input!');
    end
    if matIdx0(1) > 0
        matIdx = matIdx0(1);
    end
    
end

if matIdx ~= 0
    % determine the label of the matrix
    param.mat = obj.matrix.label{matIdx};
end

% get the point group symmetry operators of the center position
% (coupling/atom)
pOp = swsym.point(obj.lattice.sym,center(:,1));

% convert the matrices to the xyz Cartesian coordinate system
A   = obj.basisvector;
pOp = mmat(A,mmat(pOp,inv(A)));

% determine the allowed matrix elements for the first coupling/anisotropy
% convert both dr and the symmetry operators into xyz coordinate system
%[aMat, aSym] = sw_basismat(pOp,dr(:,1));
[aMat, aSym] = sw_basismat(pOp,A*dr(:,1));

aMatS = aMat(:,:,~aSym);
if bondIdx
    aMatA = aMat(:,:,aSym);
else
    aMat  = aMatS;
    aMatA = zeros(3,3,0);
    aSym  = false(1,9);
end

nSymMat = size(aMat,3);

% convert aMat in case subIdx is non-zero for coupling matrices
mod_mat = false;
if param.subIdx > 1 && anisoIdx == 0 && gIdx == 0
    mod_mat = true;
    % get the matrices for the first bond
    % TODO use cache.symop
    param0.subIdx = 1;
    param0.pref   = [];
    f0 = obj.fileid;
    obj.fileid(0);
    aMat0 = obj.getmatrix(param0);
    obj.fileid(f0);
    
    % save the original matrix of the object
    mat0 = obj.matrix.mat(:,:,matIdx);
    
    tMat = zeros(3,3,nSymMat);
    for ii = 1:nSymMat
        obj.matrix.mat(:,:,matIdx) = aMat0(:,:,ii);
        tMat(:,:,ii) = obj.couplingtable(bondIdx).matrix(:,:,param.subIdx);
    end
        
    aMat1 = zeros(3,3,nSymMat);
    for ii = 1:nSymMat
        aMat1(:,:,ii) = sum(bsxfun(@times,aMat0,permute(linsolve(reshape(tMat,9,[]),reshape(aMat(:,:,ii),[9 1])),[2 3 1])),3);
    end
    
    obj.matrix.mat(:,:,matIdx) = mat0;
    aMat = aMat1;
end

dVect = permute([aMatA(2,3,:) aMatA(3,1,:) aMatA(1,2,:)],[2 3 1]);

% Create the prefactors
if iscell(param.pref)
    pref = zeros(1,nSymMat);
    if mod(numel(param.pref{1}),2) == 0
        % create the proper prefactor vector
        pref(param.pref{1}(1:2:end)) = param.pref{1}(2:2:end);
    elseif numel(param.pref{1}) == 1
        % create Heisenberg coupling
        pref = param.pref{1};
    elseif numel(param.pref{1}) == 3
        % use prefactors for the antisymmtric matrices only
        % select antisymmetric matrices
        if sum(aSym) == 0
            error('spinw:getmatrix:NoAsym','No asymmmetric matrix is allowed by symmetry!');
        elseif sum(aSym) > 3
            error('spinw:getmatrix:SymError','Error of point group symmetry!');
        elseif sum(aSym) < 3
            warning('spinw:getmatrix:Asym','Less than 3 assymetric matrices are allowed by symmetry!');
        end
        pref(aSym) = param.pref{1}(1:sum(aSym));
    else
        error('spinw:getmatrix:WrongInput','Wrong value of pref, see help!');
    end
    param.pref = pref;
end

if nargout == 0
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
                eStr{ii} = [eStr{ii} '-D' char(48+jj)];
            elseif abs(dVect(ii,jj)-1) < tol
                % +1
                eStr{ii} = [eStr{ii} 'D' char(48+jj)];
            else
                % other number
                eStr{ii} = [eStr{ii} num2str(dVect(ii,jj)) 'D' char(48+jj)];
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
    if bondIdx
        fprintf('\nThe symmetry analysis of the coupling between atom %d and atom %d:\n',atom1(1),atom2(1));
        fprintf(' lattice translation vector: [%d,%d,%d]\n',dl(:,1));
        fprintf(' distance: %5.3f Angstrom\n',norm(obj.basisvector*dr(:,1)));
        fprintf(' center of bond (in lattice units): [%5.3f,%5.3f,%5.3f]\n', center(:,1));
    elseif anisoIdx
        fprintf('\nThe symmetry analysis of the anisotropy matrix of atom %d (''%s''):\n',anisoIdx,obj.unit_cell.label{anisoIdx});
        fprintf(' position (in lattice units): [%5.3f,%5.3f,%5.3f]\n', center(:,1));
    else
        fprintf('\nThe symmetry analysis of the g-tensor of atom %d (''%s''):\n',gIdx,obj.unit_cell.label{gIdx});
        fprintf(' position (in lattice units): [%5.3f,%5.3f,%5.3f]\n', center(:,1));        
    end
    if ~isempty(param.mat)
        fprintf(' label of the assigned matrix: ''%s''\n',param.mat);
    end
    
    fprintf(' allowed elements in the symmetric matrix:\n');
    fprintf('  S = %s\n',smatStr);
    if bondIdx
        fprintf(' allowed components of the Dzyaloshinskii-Moriya vector:\n');
        fprintf('  D = %s\n\n',amatStr);
    end
    
    if mod_mat
        fprintf('Be carefull, the output matrices are corresponding to subIdx = 1!\n');
    end
end

% calculate interaction matrix using the prefactors
if numel(param.pref) > 1
    if numel(param.pref)~=size(aMat,3)
        error('spinw:getmatrix:WrongInput','Wrong number of elements in pref option!');
    end
    param.pref = param.pref(:);
    aMat = sum(repmat(permute(param.pref,[2 3 1]),[3 3 1]).*aMat,3);
elseif numel(param.pref) == 1
    aMat = eye(3)*param.pref;
end

if nargout>0
    aMatOut  = aMat;
    param.matIdx = matIdx;
    paramOut = param;
    pOpOut  = pOp;
end

end