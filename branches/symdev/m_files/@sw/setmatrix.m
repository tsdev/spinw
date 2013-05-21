function setmatrix(obj, varargin)
% changes the selected matrix of sw object.
%
% setmatrix(obj, 'Option1', Value1, ...)
%
% Options:
%
% One of the below options has to be given:
%
% label         Label of the matrix that is already assigned to either as
%               anisotropy or coupling only once.
% mat_idx       Index of the matrix, stored in obj.matrix. Alternative to
%               the 'label' option.
% coupling_idx  Value of the obj.coupling.idx, that defines the coupling,
%               for which the symmetry allowed matrices has to be
%               determined.
% aniso_idx     Value of the obj.matom.idx, that selects a magnetic atoms,
%               for which the symmetry allowed anisotropy matrices has to
%               be determined.
%
% Optional inputs:
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
%           For anisotropy matrices antisymmetric matrices are not allowed.
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
% See also SW.GETMATRIX.
%

if nargin == 1
    help sw.setmatrix;
    return;
end

[aMat, param] = obj.getmatrix(varargin{:});

if isempty(param.pref)
    % Heisenberg coupling is always allowed by symmetry!
    aMat = eye(3);
end

if param.mat_idx > 0
    obj.matrix.mat(:,:,param.mat_idx) = aMat;
else
    error('sw:setmatrix:WrongInput','It is not possible to unambiguously select a matrix from the input options!')
end

end