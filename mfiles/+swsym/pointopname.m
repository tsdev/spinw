function [symStr, info] = pointopname(symOp,r)
% determines the name of the point operators
%
% [opStr, info] = SWSYM.POINTOPNAME(symOp,{r})
%
% Input:
%
% O         Stack of symmetry operators in a matrix with dimesions of
%           [3,4,nOp].
% r         Position vector where the point group operators are evaluated
%           (if symOp contains translation). Default value is [0 0 0].
%
% Output:
%
% opStr     The name of the operators in a cell. Only the unique names are
%           kept. The names used are:
%               'E' 'i' 'sigma' 'C2' 'C3' 'C4' 'C6' 'S3' 'S4' 'iC3'
%
% Example:
% Determines the point group operators of space group 'P -3'
%
%   O = swsym.operator('P -3');
%   swsym.pointopname(O)
%

if nargin < 2
    r = [0;0;0];
end

symOp = swsym.point(symOp,r);
%O = O(1:3,1:3,:);

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

tr0   = [3 -3  1 -1  0  1  2  0 -1 -2];
N0    = [1  2  2  2  3  4  6  6  4  6];
label = {'1' '-1' 'm' '2' '3' '4' '6' '-3' '-4' '-6'};

nOp = size(symOp,3);
tr = zeros(1,nOp);
N  = zeros(1,nOp);

for ii = 1:nOp
    % trace
    tr(ii) = sum(diag(symOp(:,:,ii)));
    % order
    N(ii) = swsym.oporder([symOp(:,:,ii) zeros(3,1)]);
end

[isSym, sIdx] = ismember([tr;N]',[tr0;N0]','rows');

if ~all(isSym)
    warning('pointopname:Error','not all symmetry is known!')
end


[symStr,~,idx] = unique(sIdx);

symStr = label(symStr);

[multi] = accumarray(idx,idx*0+1);

for ii = 1:numel(symStr)
    if multi(ii)>1
        symStr{ii} = [num2str(multi(ii)) symStr{ii}];
    end
end

if nargout>1
    info.name  = symStr;
    info.multi = multi;
    info.type  = unique(sIdx);
end

end