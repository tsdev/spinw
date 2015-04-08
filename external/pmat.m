function Mpad = pmat(M,dims,pVal)
% pads matrix to the required dimensions with a given value
%
% Mpad = PMAT(M,dims,{pVal})
%
% The output M is padded with zeros or the value given by pVal. If dims
% define a smaller dimension than the input matrix, the matrix gets cutted.
%
% Input:
%
% M         Input matrix with arbitrary dimensions.
% dims      Row vector, gives the dimensions of the output matrix.
% pVal      Value that is used for padding, default is 0.
%
% Output:
%
% Mpad      Output matrix with the padded values.
%

% create empty matrix with the padded value
Mpad = zeros([dims 1])+pVal;

% dimensions of the input matrix
dimi = size(M);
if numel(dimi)<numel(dims)
    dimi(end+1:numel(dims)) = 1;
elseif numel(dims)<numel(dimi)
    dims(end+1:numel(dimi)) = 1;
end

% common dimensions of the input matrix and the output
mindims = min(dims,dimi);

% generate the common indices in a cell
minIdx = arrayfun(@(x)1:x,mindims,'UniformOutput',false);

% copy matrix to the padded matrix
Mpad(minIdx{:}) = M(minIdx{:});

end