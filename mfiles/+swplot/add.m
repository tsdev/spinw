function add(hAdd, hFigure)
% adds a graphical object to the hgtransform of an swplot figure
%
% SWPLOT.ADD(hAdd, {hFigure})
%
% It adds a graphical object to the hgtransform object of the figure to
% enable continuous rotation with the mouse.
%
% Input:
%
% hAdd      Either vector of the handles of the graphical objects, or
%           struct with dimensions of [1 nObject] with a handle field each
%           contains a graphical object handle. The struct can contain any
%           number of the following fields as well:
%               'name'      Default value is 'general' if not given. The
%                           name identifies groups of objects.
%               'text'      Text that is shown in the tooltip when clicking
%                           on the object.
%               'position'  Position of the object, see swplot.plot for
%                           details.
%               'label'     Label that is shown in the legend.
%               'legend'    Type of legend, see swplot.plot for details.
%               'type'      Type of graphical object, see swplot.plot.
% hFigure   The handle of the figure (number in the figure title bar). The
%           default is the active swplot figure if the argument is not
%           provided by the user or it is empty matrix.
%
% See also SWPLOT, SWPLOT.FIGURE, HGTRANSFORM, SWPLOT.PLOT.
%

if nargin == 0
    help swplot.add
    return
end

if nargin < 2 || isempty(hFigure)
    % find active figure
    hFigure = swplot.activefigure;
end

hAxis = getappdata(hFigure,'axis');

cva = get(hAxis,'CameraViewAngle');
hTransform = getappdata(hFigure,'h');

% fields of struct to store objects
fNames = {'handle' 'number' 'name' 'type' 'label' 'position' 'text' 'legend'};

if isappdata(hFigure,'objects')
    sObject = getappdata(hFigure,'objects');
else
    c0 = cell(1,0);
    sInit = [fNames; repmat({c0},[1 numel(fNames)])];
    sObject = struct(sInit{:});
end

% keep the additional fields
fNames = fNames(3:end);

% find the maximum element value
if isempty(sObject)
    nMax = 0;
else
    nMax = round(max([sObject(:).number]));
end

% convert a simple list of handles to the required structure to store in
% swplot
nObjAdd = numel(hAdd);

if ~isstruct(hAdd)
    hAddC = num2cell(hAdd);
    hAdd = struct('handle',cell(1,nObjAdd));
    [hAdd.handle] = hAddC{:};
else
    hAddS = hAdd;
    if ~isfield(hAddS,'handle')
        error('add:WrongInput','handle field is required!');
    end
    hAdd = struct('handle',cell(1,nObjAdd));
    [hAdd(:).handle] = hAddS(:).handle;
    
    for ii = 1:numel(fNames)
        % copy only the allowed field names
        if isfield(hAddS,fNames{ii})
            [hAdd(:).(fNames{ii})] = hAddS(:).(fNames{ii});
        end
    end
end

% .number
num1 = num2cell(nMax+(1:nObjAdd));
[hAdd(:).number] = num1{:};

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
legend0 = [1 3 1 1 0 0 1];
type0   = {'arrow' 'ellipsoid' 'cylinder' 'circle' 'line' 'text' 'facepatch'};
% create dictionary to convert string to number
K       = containers.Map(type0,1:numel(type0));

for ii = 1:numel(fNames)
    if ~isfield(hAdd,fNames{ii})
        switch fNames{ii}
            case 'name'
                % .type
                defval1 = repmat({'general'},nObjAdd,1);
            case 'position'
                % .position
                defval1 = repmat({nan(3,2)},nObjAdd,1);
            case 'label'
                % .label, default is empty label
                defval1 = repmat({''},nObjAdd,1);
            case 'text'
                % .text
                defval1 = {hAdd(:).label};
            case 'type'
                % .type, default is taken from tag/type property
                defval1 = get([hAdd(:).handle],'Tag');
                type2   = get([hAdd(:).handle],'Type');
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
                    lIdx    = cell2mat(values(K,{hAdd(:).type}));
                    defval1 = num2cell(legend0(lIdx));
                catch
                    % .legend, default colored rectangle for all objects if
                    % above method fails
                    warning('add:WrongInput','The type of graphical object added to the plot might be not fully supported!');
                    defval1 = repmat({1},nObjAdd,1);
                end
        end
        [hAdd(:).(fNames{ii})] = defval1{:};
    end
end

% create faceindex data for facepatch type object
hNew = [hAdd(:).handle];
% find face patch objects
hPatch   = getappdata(hFigure,'facepatch');
facePatchIdx = find(hPatch==hNew);
nFacePatch   = numel(facePatchIdx);

% find edgepatch objects
hEPatch      = getappdata(hFigure,'edgepatch');
edgePatchIdx = find(hEPatch==hNew);
nVertexPatch   = numel(edgePatchIdx);

if nFacePatch > 0
    % lets check if all new hPatch handles are the same type of object
    type0 = hAdd(facePatchIdx(1)).type;
    if ~all(ismember({hAdd(facePatchIdx).type},type0))
        error('add:WrongInput','All patch objects have to be the same type!');
    else
        fIdx = getappdata(hPatch,'facenumber');
        nI = size(fIdx,1);
        F  = get(hPatch,'Faces');
        nF = size(F,1);
        
        nFacePerObject = (nF-nI)/nFacePatch;
        if ceil(nFacePerObject)-nFacePerObject > 0
            error('add:WrongInput','All patch objects have to be the same type!');
        end
        
        fIdxAdd = repmat([hAdd(facePatchIdx).number],[nFacePerObject 1]);
        fIdx = [fIdx; fIdxAdd(:)];
        setappdata(hPatch,'facenumber',fIdx);
        
    end
end

if nVertexPatch > 0
    % lets check if all new hPatch handles are the same type of object
    type0 = hAdd(edgePatchIdx(1)).type;
    if ~all(ismember({hAdd(edgePatchIdx).type},type0))
        error('add:WrongInput','All patch objects have to be the same type!');
    else
        fIdx = getappdata(hEPatch,'vertexnumber');
        nI = size(fIdx,1);
        V  = get(hEPatch,'Vertices');
        nV = size(V,1);
        
        nVertexPerObject = (nV-nI)/nVertexPatch;
        if ceil(nVertexPerObject)-nVertexPerObject > 0
            error('add:WrongInput','All patch objects have to be the same type!');
        end
        
        fIdxAdd = repmat([hAdd(edgePatchIdx).number],[nVertexPerObject 1]);
        fIdx = [fIdx; fIdxAdd(:)];
        setappdata(hEPatch,'vertexnumber',fIdx);
        
    end
end


% add all the objects to the hgtransform object except the facepatch
% handles that are already registered
hNew([facePatchIdx edgePatchIdx]) = [];
if ~isempty(hTransform)
    set(hNew,'Parent',hTransform);
end
set(hNew,'Clipping','Off');

% add callback function for showing the tooltips
%set([hNew hPatch],'ButtonDownFcn',@(obj,hit)swplot.tooltipcallback(obj,hit,hFigure,hTransform));
swplot.tooltip(swplot.tooltip,hFigure);


% comb together the handles of the old and new graphical objects.
sObject = [sObject hAdd(:)'];

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
set(hAxis,'CameraViewAngle',cva);
material('shiny');

end