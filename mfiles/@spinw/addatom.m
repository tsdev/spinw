function addatom(obj, varargin)
% adds new atom to an spinw object
%
% ADDATOM(obj,'Option1', Value1, ...)
%
% Input:
%
% obj       spinw object
%
% Options:
%
%   r       Atomic positions, dimensions are [3 nAtom]. No default value!
%   S       Spin of the atoms, dimensions are [1 nAtom], for non-magnetic
%           atoms set S to zero. Default spin is generated from the given
%           label of the atom. For example if 'label' is 'MCr3+' or 'Cr3+'
%           then the high spin of S=3/2 is automatically generated. The
%           high spin values for every ion is stored in the last column of
%           the ion.dat file. If the atom type is unknown S=0 is assumed.
%   label   Names of the atoms for plotting and form factor
%           calculations (see ion.dat), it is a cell, optional.
%           Example:
%           {'atom1' 'atom2' 'atom3'}
%           Default value is 'atomI', where I is the atom index.
%   color   Colors of the atoms for plotting, dimensions are [3 nAtom],
%           where each column describes an RGB color. Each value is between
%           0 and 255, optional. Default value is [255;165;0] for each
%           atom.
%           Alternatively a name of the color can be given as a string, for
%           example 'White', for multiple atoms package it into a cell. For
%           the list of colors, see sw_colorname().
%
% Output:
%
% The function creates extra elements in the 'unit_cell' field of the obj
% spinw object.
%
% Example:
%
% crystal.ADDATOM('r',[0 1/2; 0 0; 0 0],'S',[1 0],'color',{'red' 'blue'})
%
% Adds a magnetic atom (S=1) at position (0,0,0) and a non-magnetic one at
% (1/2 0 0) with red and blue color respectively.
%
% See also SPINW.GENLATTICE, SPINW.ADDMATRIX, SW_COLORNAME.
%

if nargin < 2
    help sw.addatom;
    return;
end

inpForm.fname  = {'r'     'S'      'label' 'color' };
inpForm.defval = {[]      []       {}      []      };
inpForm.size   = {[-1 -2] [-3 -4]  [-5 -6] [-7 -8] };
inpForm.soft   = {false   true     true    true    };

newAtom = sw_readparam(inpForm, varargin{:});

% number of old atoms
nOldAtom = numel(obj.unit_cell.S);

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

for ii = 1:numel(newAtom)
    if  ~iscell(newAtom(ii).label)
        newAtom(ii).label = {newAtom(ii).label};
    end
end

% Generate atom colors.
if ~isfield(newAtom,'color') || isempty([newAtom.color])
    for ii = 1:numel(newAtom)
        for jj = 1:numel(newAtom(ii).r)/3
            newAtom(ii).color(:,jj) = sw_atomdata(newAtom(ii).label{1},'color');
        end
    end
else
    for ii = 1:numel(newAtom)
        newAtom(ii).color = sw_colorname(newAtom(ii).color);
    end
    
end

for ii = 1:numel(newAtom)
    if ~any(size(newAtom(ii).r)-[1 3])
        newAtom(ii).r = newAtom(ii).r';
    end
    
    % Generate spins, default is generated from the label of the atom.
    if isempty(newAtom(ii).S)
        for jj = 1:numel(newAtom(ii).r)/3
            [~,~,S0] = sw_mff(newAtom(ii).label{jj});
            newAtom(ii).S(jj) = S0;
        end
    end
    
    if size(newAtom(ii).S,2) == 1
        newAtom(ii).S = newAtom(ii).S';
    end
    
    if obj.symbolic
        if ~isa(newAtom(ii).S,'sym')
            symS = sym('');
            for jj = 1:numel(newAtom(ii).r)/3
                if newAtom(ii).S(jj) > 0
                    symS(jj) = sym(['S_' num2str(jj+nOldAtom)],'positive');
                else
                    symS(jj) = sym(0);
                end
            end
            newAtom(ii).S = symS;
        end
    end
    
    if size(newAtom(ii).label,2) == 1
        newAtom(ii).label = newAtom(ii).label';
    end
    
    newAtom(ii).color = int32(newAtom(ii).color);
    newObj.unit_cell  = newAtom(ii);
    validate(newObj,'unit_cell');
    
    % check whether atom exists already
    occupied = false;
    for jj = 1:size(newObj.unit_cell.r,2)
        if any(sum(abs(bsxfunsym(@minus,obj.unit_cell.r,newObj.unit_cell.r(:,jj))),1) < 1e-2)
            occupied = true;
        end
    end
    
    if occupied
        warning('sw:addatom:AtomExists','Atomic position is already occupied, new atoms skipped!');
        return;
    end
    
    obj.unit_cell.r     = [obj.unit_cell.r     mod(newObj.unit_cell.r,1)];
    obj.unit_cell.S     = [obj.unit_cell.S     newObj.unit_cell.S];
    obj.unit_cell.label = [obj.unit_cell.label newObj.unit_cell.label];
    obj.unit_cell.color = [obj.unit_cell.color newObj.unit_cell.color];
end

validate(struct(obj));

end