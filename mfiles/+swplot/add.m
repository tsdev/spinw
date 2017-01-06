function add(hAdd, hFigure, name)
% adds a graphical object to a figure
%
% SWPLOT.ADD(hAdd, {hFigure},{name})
%
% It adds a graphical object to the hgtransform object of the figure for
% continuous rotation with the mouse.
%
% Input:
%
% hAdd      Either vector of the handles of the graphical objects, or
%           struct with dimensions of [1 noBject] with a handle field each
%           contains a graphical object handle. The struct can contain any
%           number of the following fields as well:
%               'text'
%               'position'
%               'label'
%               'legend'
%               'type'
% hFigure   The handle of the figure (number in the figure title bar). The
%           default is the active swplot figure if the argument is not
%           provided by the user or it is empty matrix.
% name      String, the name of the objects. It can be used for finding the
%           object handles after plotting.
%
% See also SWPLOT, SWPLOT.FIGURE, HGTRANSFORM.
%

if nargin == 0
    help swplot.add
    return
end

if nargin < 2 || isempty(hFigure)
    % find active figure
    hFigure = swplot.activefigure;
end

if nargin < 3
    name = 'general';
end

cva = get(gca,'CameraViewAngle');
hTransform = getappdata(hFigure,'h');

if isappdata(hFigure,'objects')
    sObject = getappdata(hFigure,'objects');
else
    sObject = struct;
end

% get the largest existing number in sObject
namesObj = fieldnames(sObject);
nMax = 0;
for ii = 1:numel(namesObj)
    nMax = max([nMax sObject.(namesObj{ii})(:).number]);
end

% convert a simple list of handles to the required structure to store in
% swplot
nObjAdd = numel(hAdd);

fNames = {'type' 'label' 'position' 'text' 'legend'};

if ~isstruct(hAdd)
    hAddC = num2cell(hAdd);
    hAdd = struct;
    hAdd.(name) = struct('handle',cell(1,nObjAdd));
    [hAdd.(name).handle] = hAddC{:};
else
    hAddS = hAdd;
    hAdd  = struct;
    if ~isfield(hAddS,'handle')
        error('add:WrongInput','handle field is required!');
    end
    hAdd.(name) = struct('handle',cell(1,nObjAdd));
    [hAdd.(name)(:).handle] = hAddS(:).handle;

    for ii = 1:numel(fNames)
        % copy only the allowed field names
        if isfield(hAddS,fNames{ii})
            [hAdd.(name)(:).(fNames{ii})] = hAddS(:).(fNames{ii});
        end
    end
end

% .number
num1 = num2cell(nMax+(1:nObjAdd));
[hAdd.(name)(:).number] = num1{:};

% define default legend type:
% 0     no legend 
% 1     colored rectangle
% 2     dashed rectangle
% 3     colored sphere
% for the following types
% 1 'arrow'
% 2 'ellipsoid'
% 3 'cylinder'
% 4 'circle'
% 5 'line'
% 6 'text'
legend0 = [1 3 1 1 0 0];
type0   = {'arrow' 'ellipsoid' 'cylinder' 'circle' 'line' 'text'};
% create dictionary to convert string to number
K       = containers.Map(type0,1:6);

for ii = 1:numel(fNames)
    if ~isfield(hAdd.(name),fNames{ii})
        switch fNames{ii}
            case 'position'
                % .position
                defval1 = repmat({nan(3,2)},nObjAdd,1);
            case 'label'
                % .label, default is empty label
                defval1 = repmat({''},nObjAdd,1);
            case 'text'
                % .text
                defval1 = {hAdd.(name)(:).label};
            case 'type'
                % .type, default is taken from tag/type property
                defval1 = get([hAdd.(name)(:).handle],'Tag');
                type2 = get([hAdd.(name)(:).handle],'Type');
                if ~iscell(defval1)
                    defval1 = {defval1};
                end
                if ~iscell(type2)
                    type2 = {type2};
                end
                tIdx  = cellfun(@(C)isempty(C),defval1);
                defval1(tIdx) = type2(tIdx);
            case 'legend'
                try
                    lIdx    = cell2mat(values(K,{hAdd.(name)(:).type}));
                    defval1 = num2cell(legend0(lIdx));
                catch
                    % .legend, default colored rectangle for all objects if
                    % above method fails
                    warning('add:WrongInput','The type of graphical object added to the plot might be not fully supported!');
                    defval1 = repmat({1},nObjAdd,1);
                end
        end
        [hAdd.(name)(:).(fNames{ii})] = defval1{:};
    end
end

% add all the objects to the hgtransform object.
hSelect = [hAdd.(name)(:).handle];
if ~isempty(hTransform)
    set(hSelect,'Parent',hTransform);
end
set(hSelect,'Clipping','Off');

% add callback function for showing the tooltips
for ii = 1:nObjAdd
    str0 = hAdd.(name)(ii).text;
    set(hAdd.(name)(ii).handle,'ButtonDownFcn',@(i,j)swplot.tooltip(str0,hFigure));
end
    
% comb together the handles of the old and new graphical objects.
if isfield(sObject,name)
    sObject.(name) = [sObject.(name) hAdd.(name)];
else
    sObject.(name) = hAdd.(name);
end

% Shift the origin to center the plot.
if isappdata(hFigure,'param')
    param = getappdata(hFigure,'param');
else
    param = struct;
end

% center object if it is crystal and hgtransform exists
% TODO change to BV matrix
if isfield(param,'range') && isappdata(hFigure,'obj') && ~isempty(hTransform)
    range       = param.range;
    basisVector = getappdata(hFigure,'obj');
    basisVector = basisVector.basisvector;
    T           = makehgtform('translate',-sum(basisVector * sum(range,2)/2,2)');
    set(hTransform,'Matrix',get(hTransform,'Matrix')*T);
end

% Saves the object handles into the figure UserData property.
setappdata(hFigure,'objects',sObject);
setappdata(hFigure,'h',hTransform);
set(gca,'CameraViewAngle',cva);
material('shiny');

end