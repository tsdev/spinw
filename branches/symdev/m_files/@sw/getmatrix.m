function [aMat, param] = getmatrix(obj, varargin)
% gives the symmetry allowed matrices for a given coupling or anisotropy
%
% aMat = GETMATRIX(obj, 'Option1', Value1, ...)
%
% Options:
%
% One of the below options has to be given:
%
% label         Label of the matrix that is already assigned to either as
%               anisotropy or coupling only once.
% mat_idx       Index of the matrix, stored in obj.matrix. Alternative to
%               the 'label' option.
% coupling_idx  Value of the obj.coupling.idx, that defines the coupling,
%               for which the symmetry allowed matrix elements have to be
%               determined.
% aniso_idx     Value of the obj.matom.idx, that selects a magnetic atom,
%               for which the symmetry allowed anisotropy matrix elements
%               have to be determined.
% g_idx         Value of the obj.matom.idx, that selects a magnetic atom,
%               for which the symmetry allowed elemtns of the g-tensor
%               have to be determined.
%
% Optional inputs:
%
% fid           For printing of the allowed matrices and point group
%               symmetries
%               0   No text output.
%               1   Output written onto the Command Window.
%               fid Output written into a text file opened with the
%                   fid = fopen(path) command.
% tol           Tolerance for printing the output for the smallest matrix
%               element.
% pref      Defines prefactors as a vector for the symmetry allowed
%           components, dimensions are [1 nSymMat]. Alternatively, if only
%           a few of the symmetry allowed matrices have non-zero
%           prefactors, use:
%               {[6 0.1 5 0.25]}
%           This means, the 6th symmetry allowed matrix have prefactor 0.1,
%           the 5th symmetry allowed matrix have prefactor 0.25. Since
%           Heisenberg isotropic couplings are always allowed, a cell with
%           a single element will create a Heisenberg coupling, example:
%               {0.1}
%           This is identical to obj.matrix.mat = eye(3)*0.1
%           For DM interactions (antisymmetric coupling matrices), use
%           three element vector in the cell:
%               {[Dx Dy Dz]}
%           In this case, these will be the prefactors of the 3
%           antisymmetric symmetry allowed matrices. In case no crystal
%           symmetry is defined, these will define directly the components
%           of the  DM interaction in the xyz coordinate system. Be
%           carefull with the sign of the DM interaction, it depends on the
%           order of the two interacting atoms! Default value is {1}.
%           For anisotropy matrices antisymmetric matrices are not allowed.
%
% Output:
%
% aMat          If no prefactors are defined, aMat contains all symmetry
%               allowed elements of the coupling/anisotropy matrix,
%               dimensions are [3 3 nSymMat]. If prefactor is defined, it
%               is a single 3x3 matrix, that is a sum of all symmetry
%               allowed elemenets multiplied by the given prefactors.
%
% See also SW.SETMATRIX.
%

inpForm.fname  = {'label' 'mat_idx' 'aniso_idx' 'coupling_idx' 'fid' 'tol' 'pref' 'g_idx'};
inpForm.defval = {zeros(1,0) 0       0          0              0     1e-5  []     0      };
inpForm.size   = {[1 -1]  [1 1]      [1 1]      [1 1]          [1 1] [1 1] [1 -2] [1 1]  };
inpForm.soft   = {false   false      false      false          false false true   false  };

param = sw_readparam(inpForm, varargin{:});
tol = param.tol;
fid = param.fid;

if nargin == 1
    help getmatrix;
    return;
end

% Check for appropriate input
inpL = [~isempty(param.label) [param.mat_idx param.aniso_idx param.coupling_idx param.g_idx]~=0];

if sum(inpL) ~= 1
    error('sw:getmatrix:WrongInput','Exactly one of the following options have to be defined: label, mat_idx, aniso_idx, g_idx or coupling_idx!');
end

if ~isempty(param.label)
    mat_idx = find(strcmp(obj.matrix.label, param.label));
    if isempty(mat_idx)
        error('sw:setmatrix:WrongInput','Matrix label cannot be found (case sensitive)!');
    elseif numel(mat_idx) > 1
        error('sw:setmatrix:WrongInput','Multiple identical matrix labels exist!');
    else
        param.mat_idx = mat_idx;
    end
end

% Identify the coupling, anisotropy or g-tensor
if param.mat_idx ~= 0
    
    % search to which coupling/anisotropy the matrix is assigned to
    [~, selCpIdx] = find(obj.coupling.mat_idx == param.mat_idx);
    [~, selAnIdx] = find(obj.single_ion.aniso == param.mat_idx);
    [~, selGIdx ] = find(obj.single_ion.g     == param.mat_idx);
    
    Cpidx = unique(obj.coupling.idx(selCpIdx));
    Anidx = unique(obj.matom.idx(selAnIdx));
    Gidx   = unique(obj.matom.idx(selGIdx));
    
    sumNum = numel(Cpidx) + numel(Anidx) + numel(Gidx);
    
    if isempty(Cpidx) && isempty(Anidx) && isempty(Gidx)
        error('sw:setmatrix:WrongInput','Matrix is not assigned to any coupling/anisotropy/g-tensor, define idx!');
    elseif sumNum > 1
        error('sw:setmatrix:WrongInput','Matrix is assigned to multiple coupling/anisotropy/g-tensor idx, define idx!');
    else
        if ~isempty(Cpidx)
            param.coupling_idx = Cpidx;
        elseif ~isempty(Anidx)
            param.aniso_idx = Anidx;
        else
            param.g_idx = Gidx;
        end
    end
end

if param.coupling_idx
    % Coupling is defined
    iSel = obj.coupling.idx == param.coupling_idx;
    if ~any(iSel)
        error('sw:getmatrix:WrongInput','The given obj.coupling.idx does not existst, check your input!');
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
    mat_idx = obj.coupling.mat_idx(:,find(iSel,1,'first'));
    if sum(mat_idx>0) == 1
        param.mat_idx = sum(mat_idx);
    end
elseif param.aniso_idx
    % Anisotropy is defined
    center = obj.unit_cell.r(:,param.aniso_idx);
    dr     = 0;
    
    mat_idx = obj.single_ion.aniso(obj.matom.idx == param.aniso_idx);
    if isempty(mat_idx)
        error('sw:getmatrix:WrongInput','The given obj.matom.idx does not exists, check your input!');
    end
    if mat_idx(1) > 0
        param.mat_idx = mat_idx(1);
    end
else
    % g-tensor is defined
    center = obj.unit_cell.r(:,param.g_idx);
    dr     = 0;
    
    mat_idx = obj.single_ion.g(obj.matom.idx == param.g_idx);
    if isempty(mat_idx)
        error('sw:getmatrix:WrongInput','The given obj.matom.idx does not exists, check your input!');
    end
    if mat_idx(1) > 0
        param.mat_idx = mat_idx(1);
    end
    
end

if param.mat_idx ~= 0
    % determine the label of the matrix
    param.label = obj.matrix.label{param.mat_idx};
end

% get the point group symmetry operators of the center position
% (coupling/atom)
pOp = sw_pointsym(obj.lattice.sym,center(:,1));

% convert the matrices to the xyz Cartesian coordinate system
A   = obj.basisvector;
pOp = mmat(A,mmat(pOp,inv(A)));

% determine the allowed matrix elements for the first coupling/anisotropy
[aMat, aSym] = sw_basismat(pOp,dr(:,1));

aMatS = aMat(:,:,~aSym);
if param.coupling_idx
    aMatA = aMat(:,:, aSym);
else
    aMat  = aMatS;
    aMatA = zeros(3,3,0);
    aSym  = false(1,9);
end

nSymMat = size(aMat,3);

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
            error('sw:setmatrix:NoAsym','No asymmmetric matrix is allowed by symmetry!');
        elseif sum(aSym) > 3
            error('sw:setmatrix:SymError','Error of point group symmetry!');
        elseif sum(aSym) < 3
            warning('sw:setmatrix:Asym','Less than 3 assymetric matrices are allowed by symmetry!');
        end
        pref(aSym) = param.pref{1}(1:sum(aSym));
    else
        error('sw:setmatrix:WrongInput','Wrong value of pref, see help!');
    end
    param.pref = pref;
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
    if param.coupling_idx
        fprintf(fid,'\nThe symmetry analysis of the coupling between atom %d and atom %d:\n',atom1(1),atom2(1));
        fprintf(fid,' lattice translation vector: [%d,%d,%d]\n',dl(:,1));
        fprintf(fid,' distance: %5.3f Angstrom\n',norm(obj.basisvector*dr(:,1)));
        fprintf(fid,' center of bond (in lattice units): [%5.3f,%5.3f,%5.3f]\n', center(:,1));
    elseif param.aniso_idx
        fprintf(fid,'\nThe symmetry analysis of the anisotropy matrix of atom %d (''%s''):\n',param.aniso_idx,obj.unit_cell.label{param.aniso_idx});
        fprintf(fid,' position (in lattice units): [%5.3f,%5.3f,%5.3f]\n', center(:,1));
    else
        fprintf(fid,'\nThe symmetry analysis of the g-tensor of atom %d (''%s''):\n',param.g_idx,obj.unit_cell.label{param.g_idx});
        fprintf(fid,' position (in lattice units): [%5.3f,%5.3f,%5.3f]\n', center(:,1));        
    end
    if ~isempty(param.label)
        fprintf(fid,' label of the assigned matrix: ''%s''\n',param.label);
    end
    
    fprintf(fid,' allowed elements in the symmetric matrix:\n');
    fprintf(fid,'  S = %s\n',smatStr);
    if param.coupling_idx
        fprintf(fid,' allowed components of the Dzyaloshinskii-Moriya vector:\n');
        fprintf(fid,'  D = %s\n\n',amatStr);
    end
end

% calculate interaction matrix using the prefactors
if numel(param.pref) > 1
    if numel(param.pref)~=size(aMat,3)
        error('sw:getmatrix:WrongInput','Wrong number of elements in pref option!');
    end
    param.pref = param.pref(:);
    aMat = sum(repmat(permute(param.pref,[2 3 1]),[3 3 1]).*aMat,3);
elseif numel(param.pref) == 1
    aMat = eye(3)*param.pref;
end

end