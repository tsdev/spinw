function opInfo = pointopname(symOp)
% determines the names of the point group operators
%
% opInfo = SWSYM.POINTOPNAME(symOp)
%
% Input:
%
% O         Stack of point group symmetry operators in a matrix with
%           dimesions of [3,3,nOp].
%
% Output:
%
% opInfo    Structure that contains the information on the operators with
%           the following fields:
%               name    Name of the point group operator in
%                       crystallographic nomenclature. Possible names:
%                       '1','-1','m','2','3','4','6','-3','-4','-6'.
%               name2   Name in molecular symmetry nomenclature. Possible
%                       names: 'E','i','sigma','C2','C3','C4','C6','S3' 
%                       'S4','iC3'.
%               axis    Row vector, that defines the axis of the symmetry
%                       operator. It is the normal to the mirror plane and
%                       the rotation axis. Zero for identity and inversion.
%
%
% Example:
% Determines the name of the point group operators of space group 'P -3'
%
%   O = swsym.operator('P -3');
%   swsym.pointopname(O)
%

% Operator label and trace and order
%  1     --> chi =  3, N = 1
% -1     --> chi = -3, N = 2
%  m, -2 --> chi =  1, N = 2
%  2     --> chi = -1, N = 2
%  3     --> chi =  0, N = 3
%  4     --> chi =  1, N = 4
%  6     --> chi =  2, N = 6
% -3     --> chi =  0, N = 6
% -4     --> chi = -1, N = 4
% -6     --> chi = -2, N = 6

tr0    = [  3   -3    1   -1    0    1    2    0   -1   -2 ];
N0     = [  1    2    2    2    3    4    6    6    4    6 ];
axis   = [  0    0   -1    1    1    1    1   -1   -1   -1 ];
label  = { '1' '-1'  'm'  '2'  '3'  '4'  '6' '-3' '-4' '-6'};
label2 = { 'E'  'i' 'sigma' 'C2' 'C3' 'C4' 'C6' 'S3' 'S4' 'iC3'};

nOp = size(symOp,3);
tr  = zeros(1,nOp);
N   = zeros(1,nOp);

for ii = 1:nOp
    % trace
    tr(ii) = sum(diag(symOp(:,:,ii)));
    % order
    N(ii) = swsym.oporder([symOp(:,:,ii) zeros(3,1)]);
end

[isSym, sIdx] = ismember([tr;N]',[tr0;N0]','rows');

% sort operators
[~,sIdx2] = sort(sIdx);

% assign labels
label  = label(sIdx);
label2 = label2(sIdx);
axis   = axis(sIdx);

if ~all(isSym)
    warning('pointopname:Error','The given matrix is not a valid point group operator!')
end

% find the relevant axis for all operator
opAx = zeros(3,nOp);

for ii = 1:nOp
    [V,D] = eig(symOp(:,:,ii));
    if axis(ii) ~= 0
        opAx(:,ii) = V(:,(abs(diag(D)-axis(ii))<10*eps));
    end
end

% create operator info structure
opInfo = struct('name',label(:),'name2',label2(:),'axis',mat2cell(opAx',ones(1,nOp),3));
opInfo = opInfo(sIdx2);

end