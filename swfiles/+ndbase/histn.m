function [Ysum, Nsum, Cbin] = histn(X, Y, varargin)
% calculates histogram of arbitrary dimensional data
%
% function [Ysum, Nsum, Cbin] = NDBASE.HISTN(X,Y,bin1,bin2,...,'option1',value1)
%
% Any data point with NaN X or Y value is omitted from the binning. Works
% for non-uniform bins and also optimised for uniform bins with significant
% speedup.
%
% Input:
%
% X         Array with size of [nPoint nDim], that represents
%           positions in the nDim dimensional space (R^nDim).
% Y         Column vector with nPoint element, that represents values at
%           the points defined in X. {X,Y} defines a scalar field in the
%           nDim dimensional space. If Y contains multiple columns, each
%           column will be binned independently and saved into Ysum in
%           separate columns, this can represent vector field at X. The
%           dimensions of the vector field is nField, if the size of Y is
%           [nPoint nField].
% binI      Row vectors that define BIN EDGES along the I-th dimension.
%           It goes from 1 to nDim. Each vector has to have at least two
%           elements. The number of bins along the I-th dimension is equal
%           to the number of elements in the binI vector minus one. Thus:
%               nBinI = numel(binI)-1
%
% Options:
%
% emptyval  Value for empty bins in Ysum and Nsum:
%               emptyval = [emptyY emptyN]
%           Default value is [NaN 0].
%
% Output:
%
% Ysum      Array with size of [nBin1 nBin2 nBin3 ... nField]. Each pixel
%           contains the sum of the elements of Y that are within the
%           given bin boundaries:
%               bin1(i1) <= X(:,i1) < bin1(i1+1) and
%               	...
%               binN(iN) <= X(:,iN) < binN(iN+1)
%           Points of X, that are outside of the bin boundaries are
%           omitted. Empty pixels will contain NaN by default.
% Nsum      Column vector with the same number of rows as Ysum, contains
%           the number of points that are contributing to the pixel. Only
%           calculated if two output is expected. Empty pixels will contain
%           the value 0 by default. To make it possible to calculate the
%           average value of each pixel, set the default value to 1, then
%           devide Ysum with Nsum element vise:
%               Yavg = bsxfun(@rdivide,Ysum,Nsum);
% Cbin      Center bin positions, if nDim=1 the values are stored in a row
%           vector. If nDim>1 the center bin vectors are packed into a
%           cell.
%
% Example:
%
% Random points in 2D.
%
% nPoint = 1e3;
% nDim   = 2;
% bin = linspace(0,1,101);
%
% X = rand(nPoint,nDim);
% Y = sin(X(:,1)*2*pi);
% [Ysum,Nsum] = ndbase.histn(X,Y,bin,bin);
% figure
% imagesc(Ysum./Nsum);
%
%
% Create points on a square lattice.
%
% [xx,yy] = ndgrid(1:0.5:10,1:0.5:10);
%
% bin = linspace(0,11,101);
%
% [Ysum,Nsum] = ndbase.histn([xx(:) yy(:)],sin(xx(:)),bin,bin);
% figure
% imagesc(Ysum./Nsum);
%
% Oversample the sine function defined on a finite point.
% xx = 0:0.5:10;
% bin = linspace(0,11,101);
%
% Ysum = ndbase.histn([xx(:)],sin(xx(:)),bin);
% figure
% plot(Ysum,'o-');
%
% See also accumarray, interp1.

%           For points in the last bin along any dimension (I-th):
%               binI(end-1) <= X(:,I) <= binI(end)
% REMOVED

if nargin == 0
    help ndbase.histn
    return
end

% default empty values
emptyval = [NaN 0];

% get the options
if nargin>3 && ischar(varargin{end-1})
    if strcmp('emptyval',varargin{end-1})
        emptyval = varargin{end};
        
        if numel(emptyval) ~= 2
            error('ndext:histn:WrongOption','The emptyval option has to contain a two element vector!')
        end
    end
    % remove the last two input variables
    varargin = varargin(1:end-2);
end


if ~ismatrix(X)
    error(['histn:WrongInput','X requires to be an (nPoint x nDim)'...
        ' array of nPoint points in R^nDim space.']);
end

% bin vectors in a cell
bin    = varargin;
% number of dimensions of the data
N      = numel(bin);
% number of field components in Y
nField = size(Y,2);
% number of bin edges
nBin   = cellfun(@(C)numel(C),bin) - 1;

if N~=size(X,2)
    error('histn:WrongInput',['The number of bin vectors doesn''t agree '...
        'with the position space dimension R^N!'])
end

if any(nBin<1)
    error('histn:WrongInput','All bin edge vector has to have at least two elements.')
end

% determine bin indices instead of coordinate values and store in X to
% spare memory and make sure it doesn't return values larger than nBin
for ii = 1:N
    % check bin for uniformity
    dBin = diff(bin{ii}(1:2));
    if any(abs(diff(bin{ii})-dBin)>eps)
        % non uniform bin
        X(:,ii) = interp1(bin{ii},1:nBin(ii)+1,X(:,ii),'previous',nBin(ii)+1);
        % remove the elements that are equal to the last bin edge
        %X(X(:,ii)==nBin(ii)+1,ii) = NaN;
    else
        % select elements that are equal to the last bin edge
        % REMOVED
        % lastIdx = X(:,ii)==bin{ii}(end);
        % uniform bin
        X(:,ii) = floor((X(:,ii) - bin{ii}(1))/dBin)+1;
        % elements in X that are equal to the last bin edge, fit into the
        % last bin
        % REMOVED
        % X(lastIdx,ii) = bin{ii}(end);
        % make NaN for outliers
        X(X(:,ii)<1 | X(:,ii)>nBin(ii),ii) = nBin(ii)+1;
        %X(X(:,ii)>nBin(ii),ii) = NaN;
    end
end

% remove points with NaN value in either the coordinate (X) or in value(Y)
% faster to sum all bad indices than remove all NaN indices
% NaN values come from the interp1 function for outliers
nanIdx = any(isnan(Y),2);
anynan = any(nanIdx);
% create an extra fake bin to sum up pixels with NaN coordinates, it is
% faster than to call accumarray with X(~nanIdx,:)
if anynan
    %X(nanIdx,2:end) = 1;
    X(nanIdx,1) = nBin(1)+1;
end

% cell that can be used to select variable number of dimension of a matrix
selDim = repmat({':'},1,N);

% preoccupy memory for Ysum to speed up the calculation
if nField > 1
    % (kind of) nice little trick to pass variable number of arguments
    % there should be a simple deal_by_element(vect) function that would do
    % it element vise, or splitting an array into multiple output
    nBinC = num2cell(nBin+1);
    Ysum  = zeros(nBinC{:},nField);
    
    if N==1
        % treat the stupid problem that size of column vectors is Nx1, give
        % this extra 1 to feed into accumarray()
        nBin  = [nBin  0];
    end

    % the same trick but now for selecting arbitrary number of dimensions
    % to use function handle such as @mean would be very slow
    for ii = 1:nField
        Ysum(selDim{:},ii) = accumarray(X,Y(:,ii),nBin+1,[],emptyval(1));
    end
else
    if N==1
        % treat the stupid problem that size of column vectors is Nx1, give
        % this extra 1 to feed into accumarray()
        nBin  = [nBin  0];
    end

    Ysum = accumarray(X,Y,nBin+1,[],emptyval(1));
end

% cut off the extra pix again usign variable number of higher dimensions if
% there is any NaN value
selDim = cell(1,N);
for ii = 1:N
    selDim{ii} = 1:nBin(ii);
end
Ysum   = Ysum(selDim{:},:);

% calculate the number of accumulated points per bin if requested
if nargout>1
    nPoint = size(Y,1);
    Nsum = accumarray(X,ones(nPoint,1),nBin+1,[],emptyval(2));
    Nsum = Nsum(selDim{:});
end

% calculate the center bin positions
if nargout>2
    if N == 1
        Cbin = (bin{1}(1:(end-1))+bin{1}(2:end))/2;
    else
        Cbin = cellfun(@(C)(C(1:(end-1))+C(2:end))/2,bin,'UniformOutput',false);
    end
end

end