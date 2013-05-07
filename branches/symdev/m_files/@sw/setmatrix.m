function setmatrix(obj, varargin)
% changes the selected matrix of sw object.
%
% setmatrix(obj, 'Option1', Value1, ...)
%
% Options:
%
% pref      Defines prefactors as a vector for the symmetry allowed
%           components, dimensions are [1 nSymMat]. Alternatively, if only
%           a few of the symmetry allowed matrices have non-zero
%           prefactors, use:
%               {[6 0.1 5 0.25]}
%           This means, the 6th symmetry allowed matrix have prefactor 0.1,
%           the 5th symmetry allowed matrix have prefactor 0.25. Since
%           Heisenberg isotropic couplings are always allowed, a cell with
%           a single element will create a Heisenberg coupling, example:
%               {0.1}
%           This is identical to obj.matrix.mat = eye(3)*0.1
%           For DM interactions (antisymmetric coupling matrices), use
%           three element vector in the cell:
%               {[Dx Dy Dz]}
%           In this case, these will be the prefactors of the 3
%           antisymmetric symmetry allowed matrices. In case no crystal
%           symmetry is defined, these will define directly the components
%           of the  DM interaction in the xyz coordinate system. Be
%           carefull with the sign of the DM interaction, it depends on the
%           order of the two interacting atoms! Default value is {1}.
% label     String, selects one of the matrices stored in sw object. Defult
%           is '', in this case, mat_idx has to be defined.
% mat_idx   Index of the matrix, stored in sw object. Alternative to the
%           label string, default is 0.
% idx       Coupling idx value, identifies the coupling, same as the index
%           used in addcoupling(). If set to 'auto', the program tries to
%           identify the couplings, to which the selected matrix is
%           assigned - only works if the selected matrix is assigned to a
%           single coupling. Default is 'auto'.
%
% Example:
%
% setmatrix(crystal,'label','J1','pref',{[6 0.235]})
% This will set 'J1' coupling to the 6th symmetry allowed matrix, with
% prefactor 0.235.
%
% setmatrix(crystal,'label','J2','pref',{1.25})
% This will set 'J2' to antiferromagnetic Heisenberg exchange, with value
% of 1.25 meV.
%

if nargin == 1
    help sw.setmatrix;
    return;
end

inpForm.fname  = {'pref' 'label' 'mat_idx' 'idx'  'tol' };
inpForm.defval = {{1}    ''      0         'auto' 1e-5  };
inpForm.size   = {[1 -1] [1 -2]  [1 1]     [1 -2] [1 1] };

param = sw_readparam(inpForm, varargin{:});

% Identify the matrix
if isempty(param.label) && (param.mat_idx<1)
    error('sw:setmatrix:WrongInput','Either label or mat_idx has to be defined!');
end
if ~isempty(param.label)
    mat_idx = find(strcmp(obj.matrix.label, param.label));
    if isempty(mat_idx)
        error('sw:setmatrix:WrongInput','Matrix label cannot be found (case sensitive)!');
    elseif numel(mat_idx) > 1
        error('sw:setmatrix:WrongInput','Multiple identical matrix labels exist!');
    else
        param.mat_idx = mat_idx;
    end
end

% Identify the coupling
if strcmpi(param.idx,'auto')
    % search for the proper coupling idx value
    coupling = obj.coupling;
    [~, selIdx] = find(coupling.mat_idx == param.mat_idx);
    idx = unique(coupling.idx(selIdx));
    if isempty(idx)
        error('sw:setmatrix:WrongInput','Matrix is not assigned to any coupling, define idx!');
    elseif numel(idx) > 1
        error('sw:setmatrix:WrongInput','Matrix is assigned to multiple coupling idx, define idx!');
    else
        param.idx = idx;
    end
end

% generate the symmetry allowed matrices
symMat = obj.getmatrix(param.idx);
nSymMat = size(symMat,3);

% Create the prefactors
if iscell(param.pref)
    pref = zeros(1,nSymMat);
    if mod(numel(param.pref{1}),2) == 0
        % create the proper prefactor vector
        pref(param.pref{1}(1:2:end)) = param.pref{1}(2:2:end);
    elseif numel(param.pref{1}) == 1
        % create Heisenberg coupling
        pref = param.pref{1};
    elseif numel(param.pref{1}) == 3
        % use prefactors for the antisymmtric matrices only
        % select antisymmetric matrices
        aSym = find(permute(sum(sum((symMat-permute(symMat,[2 1 3])).^2,1),2),[1 3 2]) > param.tol^2);
        if numel(aSym) == 0
            error('sw:setmatrix:NoAsym','No asymmmetric matrix is allowed by symmetry!');
        elseif numel(aSym) > 3
            error('sw:setmatrix:SymError','Error of point group symmetry!');
        elseif numel(aSym) < 3
            warning('sw:setmatrix:Asym','Less than 3 assymetric matrices are allowed by symmetry!');
        end
        pref(aSym) = param.pref{1}(1:numel(aSym));
    else
        error('sw:setmatrix:WrongInput','Wrong value of pref, see help!');
    end
    param.pref = pref;
end

if numel(param.pref) > 1
    obj.matrix.mat(:,:,param.mat_idx) = obj.getmatrix(param.idx,'pref',param.pref);
else
    % Heisenberg coupling (always allowed by symmetry!)
    obj.matrix.mat(:,:,param.mat_idx) = eye(3) * param.pref;
end

end