function varargout = base(varargin)
% sets the basis vectors of an swplot figure
%
% SWPLOT.BASE(BV, {hFigure})
%
% BV is a matrix with dimensions of [3 3] and contains the three basis
% vectors of the new coordinate system as column vectors.
%
% SWPLOT.BASE(obj, {hFigure})
%
% obj is a spinw object that defines the swplot coordinate system as
% lattice units.
%
% BV = SWPLOT.BASE
%
% Returns the basis vectors stored in the swplot figure.
%
% Input:
%
% BV            Either a 3x3 matrix of the new basis vectors or a spinw
%               object where the new basis vectors will be the lattice
%               units of the stored crystal.
% hFigure       Handle of the swplot figure. Default is the selected
%               figure.
%
% See also SWPLOT.PLOT.
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