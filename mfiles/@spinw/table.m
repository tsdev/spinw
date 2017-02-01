function T = table(obj,type,index,showVal)
% outputs easy to read tables of internal data
%
% T = SPINW.TABLE(obj,type,{index},{showVal})
%
% The function only works in Matlab R2013b or newer, since the output is a
% table(), not supported in earlier versions.
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
% index     Bond index, default value is 1, to show the first neighbor
%           bonds only.
% showVal   Also show the values of the single ion terms and exchange
%           values. The values shown  are the true exchange values after
%           the symmetry operations (if there is any). Default is false.
%
% Output:
%
% T         Matlab table object.
%

if verLessThan('MATLAB','R2013b')
    warning('spinw:table:Version','This function is supported only on MATLAB R2013b or newer!')
    return
end

if nargin<3
    index = 1;
end

if nargin<4
    showVal = false;
end

switch type
    case {'matom' 'atom'}
        matom = obj.unit_cell.label(obj.matom.idx)';
        
        if isempty(matom)
            T = table;
            return
        end
        
        idx   = (1:numel(obj.matom.idx))';
        S     = obj.matom.S';
        pos   = round(obj.matom.r'*1e3)/1e3;
        
        T = table(matom,idx,S,pos);
    case 'mag'
        
        nCell = prod(double(obj.mag_str.nExt));
        
        
        matom = repmat(obj.unit_cell.label(obj.matom.idx),[1 nCell])';
        idx   = repmat(1:numel(obj.matom.idx),[1 nCell])';
        absRF = sqrt(sum(real(obj.mag_str.F(:,:,1)).^2,1));
        absIF = sqrt(sum(imag(obj.mag_str.F(:,:,1)).^2,1));
        S     = max(absRF,absIF)';
        absRF(absRF==0) = 1;
        absIF(absIF==0) = 1;
        realFhat = round(bsxfun(@rdivide,real(obj.mag_str.F(:,:,1)),absRF)'*1e3)/1e3;
        imagFhat = round(bsxfun(@rdivide,imag(obj.mag_str.F(:,:,1)),absIF)'*1e3)/1e3;
        pos   = round(obj.magtable.R'*1e3)/1e3;
        num   = (1:numel(matom))';
        T = table(num,matom(:),idx,S,realFhat,imagFhat,pos);
    case {'bond' 'coupling'}
        
        [SS,~] = obj.intmatrix('extend',false,'zeroC',true);
        
        % selector
        sel = ismember(obj.coupling.idx,index);
                
        bond   = obj.coupling.idx(sel)';
        
        if isempty(bond)
            T = table;
            return
        end
        
        dl     = double(obj.coupling.dl(:,sel))';
        idx1   = obj.coupling.atom1(sel)';
        idx2   = obj.coupling.atom2(sel)';
        length = obj.basisvector*(obj.matom.r(:,idx2)+dl'-obj.matom.r(:,idx1));
        length = round(sqrt(sum(length.^2,1))'*1e3)/1e3; % Angstrom
        matom1  = obj.unit_cell.label(obj.matom.idx(idx1))';
        matom2  = obj.unit_cell.label(obj.matom.idx(idx2))';
        mLabel = [obj.matrix.label {''}];
        matIdx = obj.coupling.mat_idx;
        matIdx(matIdx==0) = numel(mLabel);
        matrix = mLabel(matIdx(:,sel))';
        
        bondSym = logical(obj.coupling.sym(:,sel)');
        
        matrix(bondSym) = cellfun(@(C)[char(187) C],matrix(bondSym),'UniformOutput',false);
        
        if numel(matrix) == 3
            matrix = matrix(:)';
        end
        
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
            
            T = table(bond,dl,length,matom1(:),idx1,matom2(:),idx2,matrix,Jx,Jy,Jz);
        else
            T = table(bond,dl,length,matom1(:),idx1,matom2(:),idx2,matrix);
        end
        
    case {'ion' 'aniso'}
        
        [~,SI] = obj.intmatrix('extend',false,'zeroC',true,'plotmode',true);
        
        matom   = obj.unit_cell.label(obj.matom.idx)';
        
        if isempty(matom)
            T = table;
            return
        end
        
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
            T = table(matom,idx,aniso,gtensor,Ax,Ay,Az,gx,gy,gz);
        else
            T = table(matom,idx,aniso,gtensor);
        end
    case {'mat' 'matrix'}
        matrix = obj.matrix.label';
        
        if isempty(matrix)
            T = table;
            return
        end
        
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
        T =table(matrix,Mx,My,Mz,type,assigned);
        
    otherwise
        error('spinw:table:WrongInput','Wrong table type string!');
end

end