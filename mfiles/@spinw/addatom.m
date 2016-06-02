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
% r         Atomic positions, dimensions are [3 nAtom]. No default value!
% S         Spin of the atoms, dimensions are [1 nAtom], for non-magnetic
%           atoms set S to zero. Default spin is generated from the given
%           label of the atom. For example if 'label' is 'MCr3+' or 'Cr3+'
%           then the high spin of S=3/2 is automatically generated. The
%           high spin values for every ion is stored in the last column of
%           the ion.dat file. If the atom type is unknown S=0 is assumed.
% label     Names of the atoms for plotting and form factor
%           calculations (see ion.dat), it is a cell, optional.
%           Example:
%           {'atom1' 'atom2' 'atom3'}
%           Default value is 'atomI', where I is the atom index.
% color     Colors of the atoms for plotting, dimensions are [3 nAtom],
%           where each column describes an RGB color. Each value is between
%           0 and 255, optional. Default value is [255;165;0] for each
%           atom.
%           Alternatively a name of the color can be given as a string, for
%           example 'White', for multiple atoms package it into a cell. For
%           the list of colors, see sw_colorname().
% ox        Oxidation number given as a double or it will be determined
%           automatically from label. Default is 0.
% occ       Occupancy, given as double. Default is 1.
% formfact  Neutron scattering form factor, given as 9 numbers, for details
%           see the help of sw_mff().
% formfactn  Neutron scattering form factor, given as 9 numbers, for details
%           see the help of sw_mff().
% formfactx X-ray scattering form factor, given as 9 numbers, for details
%           see the help of sw_cff().
% Z         Atomic number, given as integer or determined from label
%           automatically. Default is 113 (Unobtanium).
% A         Atomic mass, given as integer. Default is -1 for the natural
%           mixture of isotopes.
% bn        Neutron scattering length, given as double.
% bx        X-ray scattering length.
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
% See also SPINW.GENLATTICE, SPINW.ADDMATRIX, SW_COLORNAME, SW_MFF, SW_CFF.
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

if isempty(newAtom.formfactn)
    newAtom.formfactn = newAtom.formfact;
end
if isempty(newAtom.bn)
    newAtom.bn = newAtom.b;
end

if ~any(size(newAtom.r)-[1 3])
    newAtom.r = newAtom.r';
end

% modulo
newAtom.r = mod(newAtom.r,1);

if size(newAtom.label,2) == 1
    newAtom.label = newAtom.label';
end

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

if isempty(newAtom.occ)
    newAtom.occ = ones(1,nNewAtom);
end

if isempty(newAtom.A)
    newAtom.A = int32(-ones(1,nNewAtom));
end

if isempty(newAtom.label)
    % Generate atom labels if not given
    if ~isempty(newAtom.Z)
        newAtom.label = sw_atomdata(newAtom.Z,'name');
    else
        newAtom.label = repmat({'atom'},[1 nNewAtom]);
        newAtom.Z     = 113+zeros(1,nNewAtom); % Unobtanium
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

% Generate spins and magnetic form factor, default is generated from the label of the atom.
if size(newAtom.S,2) == 1
    newAtom.S = newAtom.S';
end

% get form factors from label
[~,newAtom.ffn,newAtom.S0] = sw_mff(newAtom.label);
[~,newAtom.ffx]            = sw_cff(newAtom.label);

newAtom.ffn = permute(newAtom.ffn,[3 2 1]);
newAtom.ffx = permute(newAtom.ffx,[3 2 1]);

% get the auto size of the magnetic moments if not given
if isempty(newAtom.S)
    newAtom.S = newAtom.S0;
end

% get scattering form factor if given
if ischar(newAtom.formfactn)
    newAtom.formfactn = {newAtom.formfactn};
end
if iscell(newAtom.formfactn)
    [~,newAtom.ffn] = sw_mff(newAtom.formfactn);
    newAtom.ffn = permute(newAtom.ffn,[3 2 1]);
elseif ~isempty(newAtom.formfactn)
    newAtom.ffn = newAtom.formfactn;
end

% x-ray scattering form factor
if ischar(newAtom.formfactx)
    newAtom.formfactx = {newAtom.formfactx};
end
if iscell(newAtom.formfactx)
    [~,newAtom.ffx] = sw_mff(newAtom.fromfactx);
    newAtom.ffx = permute(newAtom.ffx,[3 2 1]);
elseif ~isempty(newAtom.formfactx)
    newAtom.ffx = newAtom.formfactx;
end

% include 2 zeros to make both form factor the same size
newAtom.ffn = [newAtom.ffn(1,1:8,:) zeros(1,2,size(newAtom.ffn,3)) newAtom.ffn(1,9,:)];

newAtom.ff = [newAtom.ffn;newAtom.ffx];
newAtom.b  = ones(2,nNewAtom);

newAtom.Z  = int32(newAtom.Z);


if obj.symbolic
    if ~isa(newAtom.S,'sym')
        symS = sym(0);
        for jj = 1:nNewAtom
            if newAtom.S(jj) > 0
                symS(jj) = sym(['S_' num2str(jj+nOldAtom)],'positive');
            else
                symS(jj) = sym(0);
            end
        end
        newAtom.S = symS;
    end
end

newAtom.color     = int32(newAtom.color);
newObj.unit_cell  = newAtom;
validate(newObj,'unit_cell');

cField = {'r' 'S' 'label' 'color' 'ox' 'occ' 'b' 'A' 'Z'};
for ii = 1:numel(cField)
    obj.unit_cell.(cField{ii}) = [obj.unit_cell.(cField{ii}) newObj.unit_cell.(cField{ii})];
end

obj.unit_cell.ff    = cat(3,obj.unit_cell.ff,newObj.unit_cell.ff);

validate(obj);

[~,~,rIdx] = unique(obj.unit_cell.r','rows');

% check occupancy
if any(accumarray(rIdx,obj.unit_cell.occ)>1)
    warning('spinw:addatom:WrongInput','Occupancy on some site is larger than 1!')
end

end