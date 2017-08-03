function qGrid = sw_qgrid(varargin)
% creates a Q grid
%
% qGrid = SW_QGRID('option1',value1,...)
%
% The function takes the same options as the d3d constructor in Horace. The
% following options are available:
%
% u         Row vector with 3 elements, determining the first axis in rlu.
% v         Row vector with 3 elements, determining the second axis in rlu.
% w         Row vector with 3 elements, determining the third axis in rlu.
% uoffset   Column vector with 3 elements, determining the offset of origin
%           in rlu, (fourth element is accepted but discarded).
% lab       Cell array of projection axis labels with 3 elements (4th
%           element discarded).
%
% Output:
%
% qGrid     A matrix with dimensions of [3,nAx1,nAx2,nAx3], where nAxi is
%           the number of points along axis i. Along the first dimension
%           the resiprocal space points are stored.

inpForm.fname  = {'u'     'v'     'w'     'uoffset' 'lab'                          };
inpForm.defval = {[1 0 0] [0 1 0] [0 0 1] [0 0 0]   {'(h,0,0)' '(0,k,0)' '(0,0,l)'}};
inpForm.size   = {[1 3]   [1 3]   [1 3]   [1 -1]    [1 -2]                         };

inpForm.fname  = [inpForm.fname  {'bin1' 'bin2' 'bin3' 'nExt' }];
inpForm.defval = [inpForm.defval { 0      0      0     [1 1 1]}];
inpForm.size   = [inpForm.size   {[1 -3] [1 -4] [1 -5] [1 -6] }];

param = sw_readparam(inpForm, varargin{:});

if nargin == 0
    help sw_qgrid
    return
end

% size of magnetic supercell
nExt = double(param.nExt);

% create bin vectors
bin = {param.bin1 param.bin2 param.bin3};
v   = cell(1,3);

for ii = 1:3
    bin0 = bin{ii};
    switch numel(bin0)
        case 1
            v{ii} = bin0;
        case 2
            v{ii} = bin0(1):(1/nExt(ii)):bin0(2);
        case 3
            v{ii} = bin0(1):bin0(2):bin0(3);
        otherwise
            error('sw_qgrid:WrongInput','Each bin has to be a 2 or 3 element vector!')
    end
end

% create bin vector grid
vv = zeros(numel(v{1}),numel(v{2}),numel(v{3}),3);
[vv(:,:,:,1),vv(:,:,:,2),vv(:,:,:,3)] = ndgrid(v{1},v{2},v{3});

% axis matrix, each column is an axis vector
axMat = [param.u;param.v;param.w]';

% create Q matrix
qGrid = sum(bsxfun(@times,permute(axMat,[1 3 4 5 2]),permute(vv,[5 1 2 3 4])),5);

end