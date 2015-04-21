function addcoupling(obj, matrixLabel, couplingIdx, varargin)
% assigns a predefined matrix as exchange coupling on selected bonds
%
% ADDCOUPLING(obj, matrixLabel, couplingIdx, {bondIdx})
%
% Input:
%
% obj           sw class object.
% matrixLabel   Label of the matrix, or the index.
% couplingIdx   Selects the interacting atom pairs through the coupling.idx
%               number. The coupling.idx numbers are in increasing order
%               according to the distances between magnetic atoms, for
%               example all shortest interatom distances have idx=1, second
%               shortest idx=2 and so on. couplingIdx can be vector to
%               assign the matrix to multiple inequivalent magnetic atom
%               distances.
% bondIdx       Selects the indices of bonds within coupling.idx to
%               differentiate between equal length bonds. If bondIdx
%               defined, couplingIdx has to be scalar. Optional. If the
%               crystal symmetry is not P1, bondIdx is not allowed, since
%               each equivalent coupling matrix will be calculated using
%               the symmetry operators of the space group. Optional.
% atom1         If several bonds have the same length but couple different
%               type of atoms, atom1 can restrict the selected bonds to the
%               bonds where atom1 is given by the atom label or index in
%               the sw.unit_cell list. Can be a row vector of indices or
%               cell with strings. Only works if space group is 0 (P0).
%               Default is [].
% atom2         Same as atom1, but defines the second atom in the bond.
%
% TODO
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
% See also SW, SW.GENCOUPLING, SW.ADDMATRIX.
%

if ~any(obj.atom.mag)
    error('sw:addcoupling:NoMagAtom','There is no magnetic atom in the unit cell with S>0!');
end

if isnumeric(matrixLabel)
    matrixIdx = matrixLabel;
else
    matrixIdx = find(strcmp(obj.matrix.label,matrixLabel),1,'last');
end

if isempty(matrixIdx)
    error('sw:addcoupling:CouplingTypeError','Matrix label does not exist!');
end

if nargin>3
    bondIdx = varargin{1};
    if numel(couplingIdx) > 1
        warning('sw:addcoupling:CouplingSize',['couplingIdx is '...
            'non-scalar but bondIdx is defined!']);
    end
    if obj.sym
        error('sw:addcoupling:SymmetryProblem',['bondIdx is not allowed '...
            'when crystal symmetry is not P1!']);
    end
end

warn = false;
for cSelect = 1:length(couplingIdx)
    
    index = (obj.coupling.idx==couplingIdx(cSelect));
    if ~any(index)
        error('sw:addcoupling:CouplingError',['Coupling with idx=%d does '...
            'not exist, use gencoupling with larger maxDistance and '...
            'nUnitCell parameters!'],couplingIdx(cSelect));
    end
    
    index = find(index);
    % If bondIdx is defined, it selects couplings to assign to J value.
    if exist('bondIdx','var')
        index = index(bondIdx);
    end
    
    Jmod = obj.coupling.mat_idx(:,index);
    
    for ii = 1:size(Jmod,2)
        if any(Jmod(:,ii)==matrixIdx)
            warn = true;
        elseif ~all(Jmod(:,ii))
            tIndex = find(~Jmod(:,ii),1,'first');
            Jmod(tIndex,ii) = int32(matrixIdx);
        else
            error('sw:addcoupling:TooManyCoupling',['The maximum '...
                'number of allowed couplings (3) between 2 spins are reached!']);
        end
    end
    
    obj.coupling.mat_idx(:,index) = Jmod;
    
end

if warn
    warning('sw:addcoupling:CouplingIdxWarning',['Same matrix already '...
        'assigned on some coupling!']);
end

end