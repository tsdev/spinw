function sObj = findobj(varargin)
% finds object data on swplot figure
% 
% ### Syntax
% 
% `sObj = swplot.findobj(Name,Value)`
%
% `sObj = swplot.findobj(hFigure,Name,Value)`
% 
% ### Description
% 
% `sObj = swplot.findobj(Name,Value)` finds graphical objects on the active
% swplot figure hFigure which have the given property name-value pairs. The
% possible property names are:
% * `handle`    Handle of the graphical object.
% * `objID`     Unique number of the object (increasing integer numbers).
% * `name`      Name of the object, identifies groups, such as `'atom'` for
%               all atoms.
% * `label`     Label of the objects, can identify types of atoms, etc.
%               it will accept sub strings, e.g. `'Cr'` parameter would
%               match both `'Cr1 Cr3+'` and `'Cr2 Cr3+'` labels.
%
% `sObj = swplot.findobj(hFigure,Name,Value)` search for objects on the
% swplot figure identified by the `hFigure` handle.
%
% ### Output Arguments
%
% `sObj`
% : Struct that contains all the data of the found objects.
%  
% ### See Also
%
% [swplot.delete]
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
if ~all(ismember(propName,{'handle' 'number' 'name' 'label'}))
    error('findobj:WrongInput','The given property name is not searchable!')
end

% logical indexing
pIdx = zeros(numel(propName),numel(sObj));
% find the right property
for ii = 1:numel(propName)
    switch propName{ii}
        case 'label'
            % find substrings
            pIdx(ii,:) = cellfun(@(C)~isempty(C),strfind({sObj(:).(propName{ii})},propVal{ii}));
        case 'name'
            % find exact match
            pIdx(ii,:) = ismember({sObj(:).(propName{ii})},propVal{ii});
        case {'number' 'handle'}
            pIdx(ii,:) = ismember([sObj(:).(propName{ii})],propVal{ii});
    end
end

% return the requested objects
sObj = sObj(all(pIdx,1));


end