function quickham(obj,J)
% creates magnetic Hamiltonian with a single command
%
% QUICKHAM(obj, J)
%
% The function generates the bonds from the predefined crystal structure
% and assigns exchange values to bonds such as J(1) to first neighbor, J(2)
% for second neighbor etc. The command will erase all previous bond,
% anisotropy, g-tensor and matrix definitions. Even if J(idx) == 0, the
% corresponding bond and matrix will be created.
%
% Input:
%
% obj       Spinw object.
% J         Vector containing the Heisenberg exchange values. J(1) for
%           first neighbor bonds, etc.
%

fid = obj.fileid;

obj.fileid(0);

dMax = 8;
nMax = 0;
nJ   = numel(J);

idx = 1;
% generate the necessary bonds and avoid infinite loop
while nMax < nJ && idx < 12
    obj.gencoupling('maxDistance',dMax);
    dMax = dMax+8;
    % maximum bond index
    nMax = obj.coupling.idx(end);
    idx  = idx+1;
end

obj.fileid(fid);

if nMax < nJ
    warning('The necessary bond length is too long (d>100 A), not all Js will be assigned!');
    J = J(1:nMax);
end

% clear matrix definitions
obj.matrix.mat   = zeros(3,3,0);
obj.matrix.color = int32(zeros(3,0));
obj.matrix.label = cell(1,0);

nDigit = floor(log10(nJ))+1;

for ii = 1:numel(J)
    % assign non-zero matrices to bonds
    matLabel = num2str(ii,num2str(nDigit,'J%%0%dd'));
    obj.addmatrix('label',matLabel,'value',J(ii))
    obj.addcoupling('mat',matLabel,'bond',ii)
end


end