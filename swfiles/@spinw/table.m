function T = table(obj,type,index,showVal)
% outputs easy to read tables of internal data
%
% T = SPINW.TABLE(obj,type,{index},{showVal})
%
% The function returns a table in Matlab R2013b or newer while in older
% versions a struct.
%
% For the matrix labels in the list of bonds, the '>>' sign means that the
% matrix value is determined using the symmetry operations.
%
%
% Input:
%
% obj       SpinW object.
% type      String, determines the type of data to show, values:
%               'matom'     properties of magnetic atoms in the unit cell
%               'matrix'    list of matrices
%               'ion'       single ion term in the Hamiltonian
%               'bond'      properties of selected bonds
%               'mag'       magnetic structure
% index     Indexing into the type of data to show, depending on the option
%           type:
%               'bond'      indexes the bonds (1 for first neighbors,
%                           etc.), if empty all bonds will be shown.
%               'mag'       Indexes the propagation vectors, the
%                           magnetization of the selected propagation
%                           vector will be shown.
%           Default value is 1, if empty vector ([]) is given, all
%           bonds/propagation vector will be shown.
% showVal   Also show the values of the single ion terms and exchange
%           values. The values shown  are the true exchange values after
%           the symmetry operations (if there is any). Default is false.
%
% Output:
%
% T         Matlab table object.
%

if verLessThan('MATLAB','8.2')
    isTable = false;
else
    isTable = true;
end

if nargin<3
    index = 1;
end

if nargin<4
    showVal = false;
end

varName = {};
var     = {};

switch type
    case {'matom' 'atom'}
        matom = obj.unit_cell.label(obj.matom.idx)';
        
        if ~isempty(matom)
            idx   = (1:numel(obj.matom.idx))';
            S     = obj.matom.S';
            pos   = round(obj.matom.r'*1e3)/1e3;
            
            varName = {'matom','idx','S','pos'};
            var     = {  matom,  idx,  S,  pos};
        end
    case 'mag'
        if ~isempty(obj.mag_str.F)
            nCell = prod(double(obj.mag_str.nExt));
            nK0    = size(obj.mag_str.k,2);
            
            if isempty(index)
                Fdisp    = reshape(obj.mag_str.F,3,[]);
                kIdx = repmat(1:nK0,[1 obj.nmagext]);
                kIdx = kIdx(:);
            else
                Fdisp    = obj.mag_str.F(:,:,index);
                kIdx = repmat(index,[obj.nmagext 1]);
            end
            
            matom = repmat(obj.unit_cell.label(obj.matom.idx),[1 nCell])';
            idx   = repmat(1:numel(obj.matom.idx),[1 nCell])';
            absRF = sqrt(sum(real(Fdisp).^2,1));
            absIF = sqrt(sum(imag(Fdisp).^2,1));
            S     = max(absRF,absIF)';
            absRF(absRF==0) = 1;
            absIF(absIF==0) = 1;
            realFhat = round(bsxfun(@rdivide,real(Fdisp),absRF)'*1e3)/1e3;
            imagFhat = round(bsxfun(@rdivide,imag(Fdisp),absIF)'*1e3)/1e3;
            pos   = round(obj.magtable.R'*1e3)/1e3;
            num   = (1:numel(matom))';
            if nK0>1 && nargin<3
                warning('spinw:table:Multik',['The stored magnetic structure has multiple '...
                    'proppagation vectors, showing only the first, use index to select '...
                    'different propagation vector!'])
            end
            
            if isempty(index)
                kSel = obj.mag_str.k;
            else
                kSel = obj.mag_str.k(:,index);
            end
            
            kvect = round(obj.mag_str.k(:,kIdx)'*1e5)/1e5;
            
            nKdisp = size(kSel,2);
            if nKdisp > 1
                num   = repmat(num,nKdisp,1);
                matom = repmat(matom(:),nKdisp,1);
                idx   = repmat(idx,nKdisp,1);
                pos   = repmat(pos,nKdisp,1);
            end
            
            
            if any(kSel(:))
                % show imaginary values for non-zero k-vectors
                varName = {'num' 'matom' 'idx' 'S' 'realFhat' 'imagFhat' 'pos' 'kvect'};
                var     = {num    matom   idx   S   realFhat   imagFhat   pos   kvect };
            else
                % no imag
                varName = {'num' 'matom' 'idx' 'S' 'realFhat'            'pos' 'kvect'};
                var     = {num    matom   idx   S   realFhat              pos   kvect };
            end
        end
    case {'bond' 'coupling'}
        
        [SS,~] = obj.intmatrix('extend',false,'zeroC',true);
        
        % selector
        if isempty(index)
            sel = true(size(obj.coupling.idx));
        else
            sel = ismember(obj.coupling.idx,index);
        end
        
        bond   = obj.coupling.idx(sel)';
        
        if ~isempty(bond)
            % generate subidx
            subidx = ones(size(bond));
            for ii = 2:numel(bond)
                if bond(ii)==bond(ii-1)
                    subidx(ii) = subidx(ii-1)+1;
                end
            end
            
            dl     = double(obj.coupling.dl(:,sel))';
            idx1   = obj.coupling.atom1(sel)';
            idx2   = obj.coupling.atom2(sel)';
            dr     = obj.matom.r(:,idx2)+dl'-obj.matom.r(:,idx1);
            dr     = round(dr*1e3)/1e3;
            length = obj.basisvector*(obj.matom.r(:,idx2)+dl'-obj.matom.r(:,idx1));
            length = round(sqrt(sum(length.^2,1))'*1e3)/1e3; % Angstrom
            matom1  = obj.unit_cell.label(obj.matom.idx(idx1))';
            matom2  = obj.unit_cell.label(obj.matom.idx(idx2))';
            mLabel = [obj.matrix.label {''}];
            matIdx = obj.coupling.mat_idx;
            matIdx(matIdx==0) = numel(mLabel);
            matrix = mLabel(matIdx(:,sel))';
            
            if obj.symmetry
                % show >> in from of matrix labels where symmetry operators are
                % used to generate the values
                bondSym = logical(obj.coupling.sym(:,sel)');
                matrix(bondSym) = cellfun(@(C)[char(187) C],matrix(bondSym),'UniformOutput',false);
            end
            
            if numel(matrix) == 3
                matrix = matrix(:)';
            end
            
            matom1 = matom1(:);
            matom2 = matom2(:);
            
            varName = {'idx','subidx','dl','dr','length','matom1','idx1','matom2','idx2','matrix'};
            var     = { bond, subidx,   dl, dr',   length,  matom1,  idx1,  matom2,  idx2,  matrix};
            
            if showVal
                % show the values of the matrices
                value  = zeros(3,3,numel(obj.coupling.idx));
                % exchange values
                JJ = reshape(SS.all(6:14,:),3,3,[]);
                
                % find the right exchange values
                cList = [obj.coupling.dl;obj.coupling.atom1;obj.coupling.atom2];
                jList = SS.all(1:5,:);
                [~,jSel] = ismember(jList',cList','rows');
                
                for ii = 1:numel(jSel)
                    value(:,:,jSel(ii)) = value(:,:,jSel(ii))+JJ(:,:,ii);
                end
                
                value  = permute(value(:,:,sel),[3 2 1]);
                value  = round(value*1e5)/1e5;
                Jx     = value(:,:,1);
                Jy     = value(:,:,2);
                Jz     = value(:,:,3);
                varName = [varName {'Jx','Jy','Jz'}];
                var     = [var     { Jx,  Jy,  Jz }];
            end
        end
    case {'ion' 'aniso'}
        
        [~,SI] = obj.intmatrix('extend',false,'zeroC',true,'plotmode',true);
        
        matom   = obj.unit_cell.label(obj.matom.idx)';
        
        if ~isempty(matom)
            
            idx     = (1:numel(matom))';
            mLabel  = [obj.matrix.label {''}];
            aIdx    = obj.single_ion.aniso;
            aIdx(aIdx==0) = numel(mLabel);
            aniso   = mLabel(aIdx)';
            gIdx    = obj.single_ion.g;
            gIdx(gIdx==0) = numel(mLabel);
            gtensor = mLabel(gIdx)';
            
            A = round(permute(SI.aniso,[3 2 1])*1e5)/1e5;
            g = round(permute(SI.g,[3 2 1])*1e5)/1e5;
            Ax = A(:,:,1); Ay = A(:,:,2); Az = A(:,:,3);
            gx = g(:,:,1); gy = g(:,:,2); gz = g(:,:,3);
            
            if showVal
                varName = {'matom','idx','aniso','gtensor','Ax','Ay','Az','gx','gy','gz'};
                var     = {  matom  ,idx  ,aniso  ,gtensor  ,Ax  ,Ay  ,Az  ,gx  ,gy  ,gz};
            else
                varName = {'matom','idx','aniso','gtensor'};
                var     = {  matom  ,idx  ,aniso  ,gtensor};
            end
        end
    case {'mat' 'matrix'}
        matrix = obj.matrix.label';
        
        if ~isempty(matrix)
            
            M = round(permute(obj.matrix.mat,[3 2 1])*1e5)/1e5;
            Mx = M(:,:,1); My = M(:,:,2); Mz = M(:,:,3);
            
            isBond  = ismember(1:numel(matrix),obj.coupling.mat_idx(:));
            isAniso = ismember(1:numel(matrix),obj.single_ion.aniso);
            isG     = ismember(1:numel(matrix),obj.single_ion.g);
            
            assigned = repmat({'none'},1,numel(matrix));
            assigned(isBond)  = {'bond'};
            assigned(isAniso) = {'aniso'};
            assigned(isG)     = {'gtensor'};
            assigned((isBond+isAniso+isG)>1) = {'multiple'};
            assigned = assigned';
            
            typeStr = {'Heisenberg' 'anisotropic' 'antisymmetric' 'general'};
            type = typeStr(sw_mattype(obj.matrix.mat))';
            
            varName = {'matrix','Mx','My','Mz','type','assigned'};
            var     = {  matrix  ,Mx  ,My  ,Mz  ,type  ,assigned};
        end
    otherwise
        error('spinw:table:WrongInput','Wrong table type string!');
end

% generate table or struct
if isTable
    T = table(var{:},'VariableNames',varName);
else
    T = struct;
    for ii = 1:numel(var)
        T.(varName{ii}) = var{ii};
    end
end

end