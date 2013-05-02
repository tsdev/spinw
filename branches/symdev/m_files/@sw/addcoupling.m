function obj = addcoupling(obj, matrixLabel, couplingIdx, varargin)
% assign a predefined matrix to pairs of magnetic atoms as interaction
%
% ADDCOUPLING(obj, matrixLabel, couplingIdx, {bondIdx})
%
% matrixLabel   Label of the matrix, or the index.
% couplingIdx   Selects the interacting atom pairs through the coupling.idx
%               number. The coupling.idx numbers are in increasing order
%               according to the distances between magnetic atoms, for
%               example all shortest interatom distances have idx=1, second
%               shortest idx=2 and so on. couplingIdx can be vector to
%               assign the matrix to multiple inequivalent magnetic atom
%               distances.
% {bondIdx}     Selects the indices of bonds within coupling.idx to
%               differentiate between equal length bonds. If bondIdx
%               defined, couplingIdx has to be scalar. Optional. If the
%               crystal symmetry is not P1, bondIdx is not allowed, since
%               each equivalent coupling matrix will be calculated using
%               the symmetry operators of the space group.
%

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
        warning('sw:addcoupling:CouplingSize','couplingIdx is non-scalar but bondIdx is defined!');
    end
    if obj.lattice.sym > 1
        error('sw:addcoupling:SymmetryProblem','bondIdx is not allowed when crystal symmetry is not P1!');
    end
end

warn = false;
for cSelect = 1:length(couplingIdx)
    
    index = (obj.coupling.idx==couplingIdx(cSelect));
    if isempty(index)
        error('sw:addcoupling:CouplingError','Coupling with idx=%d does not exist!',couplingIdx(cSelect));
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
            error('sw:addcoupling:TooManyCoupling','The maximum number of allowed couplings (3) between 2 spins are reached!');
        end
    end
    
    obj.coupling.mat_idx(:,index) = Jmod;
    
end

if warn
    warning('sw:addcoupling:CouplingIdxWarning','Same matrix already assigned on some coupling!');
end

end