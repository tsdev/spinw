function varargout = base(BV,hFigure)
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

if nargin < 2
    % find active figure
    hFigure = swplot.activefigure;
end

if nargin == 0 || isempty(BV)
    varargout{1} = getappdata(hFigure,'base');
    return
end

if isa(BV,'spinw')
    % get lattice from spinw object
    BV = BV.basisvector;
end

setappdata(hFigure,'base',BV);

end