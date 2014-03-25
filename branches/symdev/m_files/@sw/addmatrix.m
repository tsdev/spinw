function addmatrix(obj, varargin)
% adds new matrix that can be assigned to spins in the Hamiltonian
%
% ADDMATRIX(obj, 'Option1, Value1,...)
%
% Options:
%
% value     The value of the matrix, dimensions are  [3 3 nJ], default is
%           eye(3).
% mat       Equivalent option to 'value'.
% label     Label for plotting, strings in a cell, dimensions are [1 nJ],
%           default is 'matI', where I is the index of the matrix.
% color     Color for plotting, dimensions are  [3 nJ], default is
%           [255;0;0] for all couplings.
%
% If 'value' option is scalar, a diagonal matrix is generated with the
% given scalar value in its diagonal.
%
% Example:
% ADDMATRIX(obj,'value',eye(3))
% Adds a diagonal matrix, that can describe Heisenberg interaction.
%

if nargin < 2
    help sw.addmatrix;
    return;
end

if nargin>2
    inpForm.fname  = {'value'    'mat'      'label' 'color' };
    inpForm.defval = {[]         []         {}      []      };
    inpForm.size   = {[-1 -1 -2] [-3 -3 -4] [-5 -6] [-7 -8] };
    inpForm.soft   = {true       true       true    true    };
    
    newMat = sw_readparam(inpForm, varargin{:});
    
    if ~isempty(newMat.value)
        newMat.mat = newMat.value;
    end
    
    if numel(newMat.mat) == 1
        newMat.mat = newMat.mat*eye(3);
    end
    
    if numel(newMat.color) == 3
        newMat.color = newMat.color(:);
    end
else
    newMat = varargin{1};
end

if isa(newMat,'cell')
    newObj.mat   = newMat{1};
    newObj.label = reshape(newMat{2},1,[]);
    newObj.color = int32(newMat{3});
    newMat       = newObj;
end

if isa(newMat,'struct')
    
    % Defult Heisenberg matrix.
    if isempty(newMat(1).mat)
        if ~iscell(newMat.label)
            nLabel = 1;
        else
            nLabel = size(newMat.label,2);
        end
        nJ = max(nLabel,size(newMat.color,2));
        if nJ == 0
            error('sw:addmatrix:WrongInput','Define some matrix property!');
        else
            newMat.mat = repmat(eye(3),[1 1 nJ]);
        end
    end

    if ~isfield(newMat,'label') || isempty([newMat.label])
        idx = size(obj.matrix.mat,3)+1;
        for ii = 1:numel(newMat)
            newMat(ii).label = {};
            for jj = 1:size(newMat(ii).mat,3)
                newMat(ii).label = [newMat(ii).label {['mat' num2str(idx)]}];
                idx = idx + 1;
            end
        end
    end
    
    if ~isfield(newMat,'color') || isempty([newMat.color])
        for ii = 1:numel(newMat)
            newMat(ii).color = repmat([255;0;0],1,size(newMat(ii).mat,3));
        end
    end
    
    % If the structure is a vector concatenates it.
    newJVect     = newMat;
    newMat       = struct;
    newMat.mat   = [];
    newMat.label = {};
    newMat.color = [];
    
    for ii = 1:numel(newJVect)
        if ~any(size(newJVect(ii).color)-[1 3])
            newJVect(ii).color = newJVect(ii).color';
        end
        if ~iscell(newJVect(ii).label)
            newJVect(ii).label = {newJVect(ii).label};
        end
        
        newMat.mat   = cat(3,newMat.mat,newJVect(ii).mat);
        newMat.label = [newMat.label newJVect(ii).label];
        newMat.color = [newMat.color newJVect(ii).color];
    end
    
    % Goes through the elements one-by-one.
    for ii = 1:size(newMat.mat,3)
        newJItem.mat   = newMat.mat(:,:,ii);
        newJItem.label = newMat.label(ii);
        newJItem.color = int32(newMat.color(:,ii));
        newObj.matrix  = newJItem;
        
        validate(newObj,'matrix');
        
        cIdx = find(strcmp(obj.matrix.label,newJItem.label));
        
        if numel(cIdx)>1
            error('spinw:addmatrix:LabelError','There are multiple matrices with the same label!');
            
        elseif numel(cIdx) == 1
            % Replace the coupling type with the new one.
            if obj.symb
                obj.matrix.mat(:,:,cIdx) = newJItem.mat*sym(newJItem.label);
            else
                obj.matrix.mat(:,:,cIdx) = newJItem.mat;
            end
            obj.matrix.label(cIdx)   = newJItem.label;
            obj.matrix.color(:,cIdx) = newJItem.color;
            
            
        else
            % Add the coupling type to the list.
            if obj.symb
                for jj = 1:size(newObj.matrix.mat,3)
                    obj.matrix.mat(:,:,end+1) = newObj.matrix.mat(:,:,jj)*sym(newObj.matrix.label{jj});
                end
            else
                obj.matrix.mat = cat(3,obj.matrix.mat,newObj.matrix.mat);
            end
            obj.matrix.label = [obj.matrix.label newObj.matrix.label];
            obj.matrix.color = [obj.matrix.color newObj.matrix.color];
            
        end
    end
    validate(obj);
    
else
    error('spinw:addmatrix:SecondArgumentError','Second argument wrong type!');
end

end