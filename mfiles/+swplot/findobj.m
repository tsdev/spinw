function sObj = findobj(varargin)
% finds object data on swplot figure
%
% sObj = swplot.findobj(hFigure,'PropertyName',value,...)
%
% Finds objects on hFigure that all have the given property values. The
% possible property names are:
%
%   handle      Handle of the graphical object.
%   number      Unique number of the object (increasing integer numbers).
%   name        Name of the object, identifies groups, such as 'atom' for
%               all atoms.
%
% sObj is a struct that contains all the stored data corresponding to the
% found objects.
%
% sObj = swplot.findobj('PropertyName',value,...)
%
% Finds object on the active swplot figure.
%

if mod(nargin,2) == 0
    % find active figure
    hFigure  = swplot.activefigure;
    propName = varargin(1:2:end);
    propVal  = varargin(2:2:end);
else
    hFigure  = varargin{1};
    propName = varargin(2:2:end);
    propVal  = varargin(3:2:end);
end

% get the objects
sObj = getappdata(hFigure,'objects');

if isempty(sObj)
    return
end

% find the given properties
if ~all(ismember(propName,{'handle' 'number' 'name'}))
    error('findobj:WrongInput','The given property name is not searchable!')
end

% logical indexing
pIdx = zeros(numel(propName),numel(sObj));
% find the right property
for ii = 1:numel(propName)
    if strcmp(propName{ii},'name')
        pIdx(ii,:) = ismember({sObj(:).(propName{ii})},propVal{ii});
    else
        pIdx(ii,:) = ismember([sObj(:).(propName{ii})],propVal{ii});
    end
end

% return the requested objects
sObj = sObj(all(pIdx,1));


end