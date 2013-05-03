function addatom(obj, varargin)
% adds new atom to an sw object
%
% ADDATOM(obj,'Option1',Value1,...)
%
% Options:
%
%   r       Atomic positions, dimensions are [3 nAtom].
%   S       Spin of the atoms, dimensions are [1 nAtom], for non-magnetic
%           atoms set S to zero.
%   {label} Names of the atoms for plotting and form factor
%           calculations (see formfactor.dat), it is a cell, optional.
%           Example:
%           {'atom1' 'atom2' 'atom3'}
%   {color} Colors of the atoms for plotting, dimensions are [3 nAtom],
%           where each column describes an RGB color. Each value is between
%           0 and 255, optional.
%
% obj = ADDATOM(obj,atom)
%
% Atom is a cell, it has 4 elements, {r S label color}, where the members
% of the cell are matrices, as above.
%
% Example:
% ADDATOM(obj,'r',[0 1/2; 0 0; 0 0],'S',[1 0])
% adds a magnetic atom at position (0,0,0) and a non-magnetic one at
% (1/2 0 0).
%

if nargin < 2
    error('sw:addatom:WrongNumberOfInput','Wrong number of input!');
end

if nargin>2
    inpForm.fname  = {'r'     'S'      'label' 'color' };
    inpForm.defval = {[]      []       {}      []      };
    inpForm.size   = {[-1 -2] [-3 -4]  [-5 -6] [-7 -8] };
    inpForm.soft   = {false   false    true    true    };
    
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

if ~isfield(newAtom,'label') || isempty([newAtom.label])
    idx = 1;
    for ii = 1:numel(newAtom)
        newAtom(ii).label = {};
        for jj = 1:numel(newAtom(ii).r)/3
            newAtom(ii).label = [newAtom(ii).label {['atom' num2str(idx)]}];
            idx = idx + 1;
        end
    end
end

if ~isfield(newAtom,'color') || isempty([newAtom.color])
    for ii = 1:numel(newAtom)
        newAtom(ii).color = repmat([255;165;0],1,numel(newAtom(ii).r)/3);
    end
end

if isa(newAtom,'struct')
    for ii = 1:numel(newAtom)
        if ~any(size(newAtom(ii).r)-[1 3])
            newAtom(ii).r = newAtom(ii).r';
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