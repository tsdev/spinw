function addcoupling(obj, varargin)
% assigns a predefined matrix as exchange coupling on selected bonds
%
% ADDCOUPLING(obj, 'option1', value1, ...)
%
% Input:
%
% obj           spinw class object.
%
% Options:
%
% mat           Label or index of the matrix that will be assigned to
%               selected bonds.
% bond          Selects the interacting atom pairs through the
%               obj.coupling.idx number. The coupling.idx numbers are in
%               increasing order according to the distances between
%               magnetic atoms, for example all shortest interatom
%               distances have idx=1, second shortest idx=2 and so on.
%               bondIdx can be a vector to assign the matrix to multiple
%               inequivalent bonds.
% atom          Contains labels of atoms or index of atoms in the
%               obj.unit_cell list of atoms. If a single string label or
%               number is given, only bonds between the selected atoms will
%               be assigned. If a cell with 2 strings are given, only bonds
%               between the two selected atoms will be assigned. Only works
%               if space group is P0. Default is [].
% subIdx        If the above options are not enough to select the desired
%               bonds, using subIdx bonds can be selected one-by-one from
%               the list of bonds that fulfill the above options. Only
%               works if the space group is P0.
% type          Type of the coupling. Possible values:
%                   0       Quadratic exchange, default.
%                   1       Biquadratic exchange.
%               Can be also one of the following string: 'quadratic',
%               'biquadratic'.
%
% Output:
%
% The function adds extra entries in the 'coupling.matrix' field of the obj
% sw object.
%
% Example:
%
% ...
% cryst.addmatrix('label','J1','value',0.123)
% cryst.gencoupling
% cryst.addcoupling('J1',1)
%
% This will add the 'J1' diagonal matrix to all second shortes bonds
% between magnetic atoms.
%
% See also SPINW, SPINW.GENCOUPLING, SPINW.ADDMATRIX.
%

if ~any(obj.atom.mag)
    error('sw:addcoupling:NoMagAtom','There is no magnetic atom in the unit cell with S>0!');
end

inpForm.fname  = {'mat'  'bond' 'atom' 'subIdx' 'type'};
inpForm.defval = {[]     []     []      []      []    };
inpForm.size   = {[1 -1] [1 -2] [-3 -4] [1 -5]  [1 -6]};
inpForm.soft   = {false  false  true    true    true  };

param = sw_readparam(inpForm, varargin{:});

if ~isnumeric(param.mat)
    param.mat = find(ismember(obj.matrix.label,param.mat));
end

if isempty(param.type)
    param.type = 0*param.mat;
end

if ischar(param.type)
    param.type = {param.type};
end

if iscell(param.type)
    % check string input for the coupling type
    qSel  = strcmp(param.type,'quadratic');
    bqSel = strcmp(param.type,'biquadratic');
    param.type = nan(size(param.type));
    param.type(qSel)  = 0;
    param.type(bqSel) = 1;
    if any(isnan(param.type))
        error(['sw:addcoupling:WrongInput','Wrong coupling type, '...
            'currently only quadratic and biquadratic exchanges '...
            'are supported!']);
    end
    
end

if any(size(param.mat)~=size(param.type))
    error('sw:addcoupling:WrongInput',['A coupling type has to be '...
        'provided for each input matrix!'])
end

% select atoms
if ~isnumeric(param.atom)
    if ~iscell(param.atom)
        param.atom = {param.atom};
    end
    if numel(param.atom)>2
        error('sw:addcoupling:WrongInput','Only two different atom label can be given at once!');
    end
    aIdx1 = find(ismember(obj.unit_cell.label,param.atom{1}));
    if numel(param.atom)>1
        aIdx2 = find(ismember(obj.unit_cell.label,param.atom{2}));
    else
        aIdx2 = aIdx1;
    end
end

if ~isempty(param.atom) || ~isempty(param.subIdx)
    if numel(param.bond) > 1
        warning('sw:addcoupling:CouplingSize',['bond parameter has to be '...
            'scalar, only the first bond is selected!']);
        param.bond = param.bond(1);
    end
    if obj.sym
        error('sw:addcoupling:SymmetryProblem',['atom and subIdx options are not allowed '...
            'when crystal symmetry is not P0!']);
    end
end

idx = ismember(obj.coupling.idx,param.bond);
if ~any(idx)
    error('sw:addcoupling:CouplingError',['Coupling with idx=%d does '...
        'not exist, use gencoupling with larger maxDistance and '...
        'nUnitCell parameters!'],param.bond(1));
end

% Select bonds with given atoms
% convert atom indices from the unit_cell into matom indices
if ~isempty(param.atom)
    
    matom = obj.matom;
    maIdx1a = ismember(obj.coupling.atom1,find(ismember(matom.idx,aIdx1)));
    maIdx2a = ismember(obj.coupling.atom2,find(ismember(matom.idx,aIdx2)));
    maIdx1b = ismember(obj.coupling.atom1,find(ismember(matom.idx,aIdx2)));
    maIdx2b = ismember(obj.coupling.atom2,find(ismember(matom.idx,aIdx1)));
    idx = idx & ((maIdx1a & maIdx2a) | (maIdx1b & maIdx2b));
end

idx = find(idx);

% IfsubIdx is defined, subselect bonds.
if ~isempty(param.subIdx)
    idx = idx(param.subIdx);
end

if isempty(idx)
    warning('sw:addcoupling:NoBond','No matrix assigned, since no bond fulfilled the given options!')
    return
end

Jmod = obj.coupling.mat_idx(:,idx);
Tmod = obj.coupling.type(:,idx);

param.mat = int32(param.mat);
param.type   = int32(param.type);

if any(ismember(Jmod(:),param.mat))
    warning('sw:addcoupling:CouplingIdxWarning',['Same matrix already '...
        'assigned on some coupling!']);
end

if any(Jmod(3,:))
    error('sw:addcoupling:TooManyCoupling',['The maximum '...
        'number of allowed couplings (3) per bond is reached!']);
end

for ii = 1:numel(param.mat)
    idxSel = sub2ind(size(Jmod),sum(Jmod>0,1)+1,1:size(Jmod,2));
    Jmod(idxSel) = param.mat(ii);
    Tmod(idxSel) = param.type(ii);
end

obj.coupling.mat_idx(:,idx) = Jmod;
obj.coupling.type(:,idx)    = Tmod;

end