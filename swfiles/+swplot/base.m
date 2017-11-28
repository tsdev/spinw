function varargout = base(varargin)
% sets the basis vectors of an swplot figure
% 
% ### Syntax
% 
% `swplot.base(BV)`
%
% `swplot.base(obj)`
% 
% `BV = swplot.base`
%
% ### Description
% 
% `swplot.base(BV)` sets the basis vector for an swplot figure. The basis
% vectors can be used to define a non-orthogonal coordinate system for
% graphic objects.
%
% `swplot.base(obj)` sets the basis vectors to the lattice units of a given
% [spinw] object `obj`.
%  
% `BV = swplot.base` returns the basis vectors stored in the swplot figure.
%  
% 
% ### Input Arguments
% 
% `BV`
% : Either a $[3\times 3]$ matrix of the new basis vectors or a [spinw]
%   object where the new basis vectors will be the lattice
%   units and the basis vectors are generated via [spinw.basisvector].
% 
% `hFigure`
% : Handle of the [swplot] figure. Default is the active
%   figure.
% 
% ### See Also
% 
% [swplot.plot]
%

BV      = [];
hFigure = [];

switch nargin
    case 0
    case 1
        if numel(varargin{1})~=1 || isa(varargin{1},'spinw')
            BV = varargin{1};
        else
            hFigure = varargin{1};
        end
    case 2
        if numel(varargin{2})~=1 || isa(varargin{2},'spinw')
            varargin = varargin([2 1]);
        end
        BV = varargin{1};
        hFigure = varargin{2};
    otherwise
        error('base:WrongInput','Too many input arguments!')
end

if isempty(hFigure)
    % find active figure
    hFigure = swplot.activefigure;
end

if isa(BV,'spinw')
    % get lattice from spinw object
    BV = BV.basisvector;
end

if ~isempty(BV)
    setappdata(hFigure,'base',BV);
end

if isempty(BV) || nargout > 0
    varargout{1} = getappdata(hFigure,'base');
end

end