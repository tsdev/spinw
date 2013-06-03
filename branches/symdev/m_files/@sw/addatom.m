function addatom(obj, varargin)
% adds new atom to an sw object
%
% ADDATOM(obj,'Option1', Value1, ...)
%
% Options:
%
%   r       Atomic positions, dimensions are [3 nAtom]. No default value!
%   S       Spin of the atoms, dimensions are [1 nAtom], for non-magnetic
%           atoms set S to zero. Default is 1.
%   label   Names of the atoms for plotting and form factor
%           calculations (see ion.dat), it is a cell, optional.
%           Example:
%           {'atom1' 'atom2' 'atom3'}
%           Default value is 'atomI', where I is the atom index.
%   color   Colors of the atoms for plotting, dimensions are [3 nAtom],
%           where each column describes an RGB color. Each value is between
%           0 and 255, optional. Default value is [255;165;0] for each
%           atom.
%
% ADDATOM(obj,atom)
%
% Input:
%
% atom is a cell, it has 4 elements, {r S label color}, where the members
% of the cell are matrices, as above.
%
%
% Example:
%
% addatom(crystal,'r',[0 1/2; 0 0; 0 0],'S',[1 0])
% Adds a magnetic atom (S=1) at position (0,0,0) and a non-magnetic one at
% (1/2 0 0).
%

if nargin < 2
    help sw.addatom;
    return;
end

if nargin>2
    inpForm.fname  = {'r'     'S'      'label' 'color' };
    inpForm.defval = {[]      []       {}      []      };
    inpForm.size   = {[-1 -2] [-3 -4]  [-5 -6] [-7 -8] };
    inpForm.soft   = {false   true     true    true    };
    
    newAtom = sw_readparam(inpForm, varargin{:});
else
    newAtom = varargin{1};
end

if isa(newAtom,'cell')
    newObj.r     = newAtom{1};
    newObj.S     = reshape(newAtom{2},1,[]);
    newObj.label = reshape(newAtom{3},1,[]);
    newObj.color = int32(newAtom{4});
    newAtom = newObj;
end

if isa(newAtom,'struct')
    
    % Generate atom labels.
    if ~isfield(newAtom,'label') || isempty([newAtom.label])
        idx = size(obj.unit_cell.r,2)+1;
        for ii = 1:numel(newAtom)
            newAtom(ii).label = {};
            for jj = 1:numel(newAtom(ii).r)/3
                newAtom(ii).label = [newAtom(ii).label {['atom' num2str(idx)]}];
                idx = idx + 1;
            end
        end
    end
    
    % Generate atom colors.
    if ~isfield(newAtom,'color') || isempty([newAtom.color])
        for ii = 1:numel(newAtom)
            for jj = 1:numel(newAtom(ii).r)/3
                %newAtom(ii).color(jj) = repmat([255;165;0],1,numel(newAtom(ii).r)/3);
                newAtom(ii).color(:,jj) = sw_atomdata(newAtom(ii).label,'color');
            end
        end
    end
    
    for ii = 1:numel(newAtom)
        if  ~iscell(newAtom(ii).label)
            newAtom(ii).label = {newAtom(ii).label};
        end
        
        if ~any(size(newAtom(ii).r)-[1 3])
            newAtom(ii).r = newAtom(ii).r';
        end
        
        % Generate spins, default is 1.
        if isempty(newAtom(ii).S)
            newAtom(ii).S = ones(1,size(newAtom.r,2));
        end
        
        if ~any(size(newAtom(ii).color)-[1 3])
            newAtom(ii).color = newAtom(ii).color';
        end
        if size(newAtom(ii).S,2) == 1
            newAtom(ii).S = newAtom(ii).S';
        end
        if size(newAtom(ii).label,2) == 1
            newAtom(ii).label = newAtom(ii).label';
        end
        
        newAtom(ii).color    = int32(newAtom(ii).color);
        newObj.unit_cell = newAtom(ii);
        validate(newObj,'unit_cell');
        
        % check whether atom exists already
        occupied = false;
        for jj = 1:size(newObj.unit_cell.r,2)
            if any(sum(abs(bsxfun(@minus,obj.unit_cell.r,newObj.unit_cell.r(:,jj))),1) < 1e-2)
                occupied = true;
            end
        end
        
        if occupied
            warning('sw:addatom:AtomExists','Atomic position is already occupied, new atoms skipped!');
            return;
        end
        
        obj.unit_cell.r     = [obj.unit_cell.r     newObj.unit_cell.r];
        obj.unit_cell.S     = [obj.unit_cell.S     newObj.unit_cell.S];
        obj.unit_cell.label = [obj.unit_cell.label newObj.unit_cell.label];
        obj.unit_cell.color = [obj.unit_cell.color newObj.unit_cell.color];
    end
    validate(struct(obj));
    
else
    error('sw:addatom:SecondArgumentError','Second argument wrong type!');
end

end