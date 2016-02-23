function addmatrix(obj, varargin)
% adds new matrix that can be assigned to spins in the Hamiltonian
%
% ADDMATRIX(obj, 'Option1, Value1,...)
%
% Input:
%
% obj       spinw class object
%
% Options:
%
% value     The value of the matrix, dimensions are  [3 3 nJ], default is
%           eye(3). If the given value is scalar, a diagonal matrix is
%           generated with the given value in its diagonal. If the given
%           value is a 3 element vector, a DM interaction matrix is created
%           according to the following rule:
%           DM = [0 Dz -Dy;-Dz 0 Dx;Dy -Dx 0].
% mat       Alternative option name to 'value'.
% label     Label for plotting, strings in a cell, dimensions are [1 nJ],
%           default is 'matI', where I is the index of the matrix.
% color     Color for plotting, either a matrix with dimensions are  [3 nJ]
%           that contains color RGB codes (0-255), or string with the name
%           of the color (for multiple matrix the string have to be packed
%           into a cell. Default color is red.
%
% Output:
%
% The obj output will contain the added matrix in the obj.matrix field.
%
% Example:
%
% crystal.ADDMATRIX('value',eye(3))
%
% Adds a diagonal matrix, that can describe Heisenberg interaction.
%
% See also SPINW, SW_COLORNAME.
%

if nargin < 2
    help spinw.addmatrix;
    return;
end

inpForm.fname  = {'value'    'mat'      'label' 'color' };
inpForm.defval = {[]         []         {}      []      };
inpForm.size   = {[-1 -2 -3] [-4 -5 -6] [-7 -8] [-9 -10] };
inpForm.soft   = {true       true       true    true    };

newMat = sw_readparam(inpForm, varargin{:});

if ~isempty(newMat.value)
    newMat.mat = newMat.value;
end

if numel(newMat.mat) == 1
    newMat.mat = newMat.mat*eye(3);
end

if numel(newMat.mat) == 3
    M = newMat.mat;
    newMat.mat = [0 M(3) -M(2);-M(3) 0 M(1);M(2) -M(1) 0];
end


% Defult Heisenberg matrix.
if isempty(newMat(1).mat)
    if ~iscell(newMat.label)
        nLabel = 1;
    else
        nLabel = size(newMat.label,2);
    end
    nJ = nLabel;
    if nJ == 0
        error('sw:addmatrix:WrongInput','Define some matrix property!');
    else
        newMat.mat = repmat(eye(3),[1 1 nJ]);
        warning('sw:addmatrix:NoValue','No valid value was given for the new matrix, default value used!');
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
else
    for ii = 1:numel(newMat)
        newMat(ii).color = sw_colorname(newMat(ii).color);
    end
    
end

% If the structure is a vector concatenates it.
newJVect     = newMat;
newMat       = struct;
newMat.mat   = [];
newMat.label = {};
newMat.color = [];

for ii = 1:numel(newJVect)
    %if ~any(size(newJVect(ii).color)-[1 3])
    %    newJVect(ii).color = newJVect(ii).color';
    %end
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
        if obj.symb && ~isa(newJItem.mat,'sym')
            obj.matrix.mat(:,:,cIdx) = newJItem.mat*sym(newJItem.label,'real');
        else
            obj.matrix.mat(:,:,cIdx) = newJItem.mat;
        end
        obj.matrix.label(cIdx)   = newJItem.label;
        obj.matrix.color(:,cIdx) = newJItem.color;
        
    else
        % Add the coupling type to the list.
        if (obj.symb) && ~isa(newObj.matrix.mat,'sym')
            for jj = 1:size(newObj.matrix.mat,3)
                if newObj.matrix.label{jj}(end) == '-'
                    symVar = sym(newObj.matrix.label{jj}(1:(end-1)),'real');
                else
                    symVar = sym(newObj.matrix.label{jj},'real');
                end
                obj.matrix.mat = cat(3,obj.matrix.mat,newObj.matrix.mat(:,:,jj)*symVar);
            end
        else
            obj.matrix.mat = cat(3,obj.matrix.mat,newObj.matrix.mat);
        end
        obj.matrix.label = [obj.matrix.label newObj.matrix.label];
        obj.matrix.color = [obj.matrix.color newObj.matrix.color];
        
    end
end

validate(obj);

end