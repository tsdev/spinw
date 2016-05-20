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
    help spinw.addatom
    return
end

inpForm.fname  = {'r'     'S'      'label' 'color' 'ox'   'occ'  'bn'    };
inpForm.defval = {[]      []       {}      []      []     []     []      };
inpForm.size   = {[-1 -2] [-3 -4]  [1 -5] [1 -6]   [1 -7] [1 -7] [1 -7]  };
inpForm.soft   = {false   true     true    true    true   true   true    };

inpForm.fname  = [inpForm.fname  {'bx'   'formfactn' 'formfactx' 'b'    'formfact' 'A'    'Z'   }];
inpForm.defval = [inpForm.defval {[]     []          []          []     []         []     []    }];
inpForm.size   = [inpForm.size   {[1 -7] [-8 -9]     [-8 -9]     [1 -7] [-8 -9]    [1 -7] [1 -7]}];
inpForm.soft   = [inpForm.soft   {true    true       true         true  true       true   true  }];

newAtom = sw_readparam(inpForm, varargin{:});

% number of old atoms
nOldAtom = size(obj.unit_cell.r,2);
% number of new atoms
nNewAtom = numel(newAtom.r)/3;

% check oxidation number
if ~isempty(newAtom.ox)
    if numel(newAtom.ox) ~= nNewAtom
        error('spinw:addatom:WrongInput','Wrong input for option ''ox''!')
    end
end

% atomic number
if ~isempty(newAtom.Z)
    if numel(newAtom.Z) ~= nNewAtom
        error('spinw:addatom:WrongInput','Wrong input for option ''Z''!')
    end
end

% atomic mass number
if ~isempty(newAtom.A)
    if numel(newAtom.A) ~= nNewAtom
        error('spinw:addatom:WrongInput','Wrong input for option ''A''!')
    end
end

if  ~iscell(newAtom.label)
    newAtom.label = {newAtom.label};
end

if isempty(newAtom.label)
    % Generate atom labels if not given
    if ~isempty(newAtom.Z)
        newAtom.label = sw_atomdata(newAtom.Z,'name');
    else
        newAtom.label = repmat({'atom'},[1 nNewAtom]);
    end
    
    if isempty(newAtom.ox)
        newAtom.ox = zeros(1,nNewAtom);
    end
    for jj = 1:nNewAtom
        if newAtom.ox(jj)<0
            newAtom.label{jj} = [newAtom.label{jj} num2str(-newAtom.ox(jj)) '-_' num2str(jj+nOldAtom)];
        elseif newAtom.ox(jj)>0
            newAtom.label{jj} = [newAtom.label{jj} num2str(newAtom.ox(jj)) '+_' num2str(jj+nOldAtom)];
        else
            newAtom.label{jj} = [newAtom.label{jj} '_' num2str(jj+nOldAtom)];
        end
    end
elseif isempty(newAtom.Z) || isempty(newAtom.ox)
    % try to determine Z and ox values
    aOx  = zeros(1,nNewAtom);
    newAtom.name = cell(1,nNewAtom);
    
    for ii = 1:nNewAtom
        % remove the first word (good for typical .cif file labels)
        aName = strword(newAtom.label{ii},2,true);
        aName = aName{1};
        % get oxidation number
        aOx(ii) = str2double(aName(ismember(aName,'0':'9')));
        % sign
        if aName(end) == '-'
            aOx(ii) = -aOx(ii);
        end
        
        % remove numbers
        aName  = aName(ismember(aName,['A':'Z' 'a':'z']));
        % keep only the last 2 character
        if numel(aName)>2
            aName = aName(end+[-1 0]);
        end
        
        newAtom.name{ii} = aName;
    end
    
    aOx(isnan(aOx)) = 0;
    
    if isempty(newAtom.ox)
        newAtom.ox = aOx;
    end
    if isempty(newAtom.Z)
        newAtom.Z = sw_atomdata(newAtom.name,'Z');
    end
end

% Generate atom colors
if isempty(newAtom.color)
    newAtom.color = sw_atomdata(newAtom.Z,'color');
else
    newAtom.color = sw_colorname(newAtom.color);
end


if ~any(size(newAtom.r)-[1 3])
    newAtom.r = newAtom.r';
end

% Generate spins, default is generated from the label of the atom.
if isempty(newAtom.S)
    for jj = 1:nNewAtom
        [~,~,S0] = sw_mff(newAtom.label{jj});
        newAtom.S(jj) = S0;
    end
end

if size(newAtom.S,2) == 1
    newAtom.S = newAtom.S';
end

if obj.symbolic
    if ~isa(newAtom.S,'sym')
        symS = sym('');
        for jj = 1:numel(newAtom.r)/3
            if newAtom.S(jj) > 0
                symS(jj) = sym(['S_' num2str(jj+nOldAtom)],'positive');
            else
                symS(jj) = sym(0);
            end
        end
        newAtom.S = symS;
    end
end

if size(newAtom.label,2) == 1
    newAtom.label = newAtom.label';
end

newAtom.color = int32(newAtom.color);
newObj.unit_cell  = newAtom;
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

validate(struct(obj));

end