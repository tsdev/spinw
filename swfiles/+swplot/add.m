function add(hAdd, hFigure, showtooltip)
% adds a graphical object to an swplot figure
% 
% ### Syntax
% 
% `swplot.add(hAdd)`
%
% `swplot.add(hAdd,hFigure)`
% 
% ### Description
% 
% `swplot.add(hAdd)` adds a graphical object to the active swplot figure to
% enable continuous rotation with the mouse. The function adds the
% graphical objects to as a children to the [matlab.hgtransform].
%  
% `swplot.add(hAdd,hFigure)` adds the graphical objects to the figure of
% the figure handle `hFigure`.
% 
% ### Input Arguments
% 
% `hAdd`
% : Either vector of the handles of the graphical objects, or
%   struct with $n_{obj}$ number of elements with a `handle` field each
%   containing a graphical object handle. The struct can contain any subset
%   of the following fields as well:
%   * `name`      Default value is `'general'` if not given. The
%                 name identifies groups of objects.
%   * `text`      Text that is shown in the tooltip when clicking
%                 on the object.
%   * `position`  Position of the object, see [swplot.plot] for
%                 details.
%    * `label`    Label that is shown in the legend.
%    * `legend`   Type of legend, see [swplot.legend] for details.
%    * `type`     Type of graphical object, see [swplot.plot].
%    * `data`     Arbitrary data assigned to the object.
% 
% `hFigure`
% : The handle of the figure or number in the figure title. The
%   default value is the active swplot figure if `hFigure` is not given or
%   empty matrix.
% 
% ### See Also
% 
% [swplot] \| [swplot.figure] \| [matlab.hgtransform] \| [swplot.delete]
%

if nargin == 0
    swhelp swplot.add
    return
end

if nargin < 2 || isempty(hFigure)
    % find active figure
    hFigure = swplot.activefigure;
end

if nargin<3
    showtooltip = true;
end

hAxis = getappdata(hFigure,'axis');

hTransform = getappdata(hFigure,'h');

sObject = getappdata(hFigure,'objects');

% % fields of struct to store objects
fNames = {'handle' 'number' 'name' 'type' 'label' 'position' 'text' 'legend' 'data'};
%
% if isappdata(hFigure,'objects')
%     sObject = getappdata(hFigure,'objects');
% else
%     c0 = cell(1,0);
%     sInit = [fNames; repmat({c0},[1 numel(fNames)])];
%     sObject = struct(sInit{:});
% end

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
            case 'data'
                % empty .data
                defval1 = cell(1,nObjAdd);
        end
        [hAdd(:).(fNames{ii})] = defval1{:};
    end
end

% all handles
hNew = [hAdd(:).handle];
% check if single handle has multiple objects
if nObjAdd>1 && hNew(1) == hNew(2)
    % multiple objects within one patch
    
    % single handle
    hPatch = hAdd(1).handle;
    
    % lets check if all new hPatch handles are the same type of object
    type0 = hAdd(1).type;
    if ~all(ismember({hAdd(:).type},type0))
        error('add:WrongInput','All patch objects have to be the same type!');
    end
    
    if strcmp(get(hPatch,'FaceColor'),'flat')
        % face patch
        % lets check if all new hPatch handles are the same type of object
        fIdx = getappdata(hPatch,'facenumber');
        nI   = size(fIdx,1);
        F    = get(hPatch,'Faces');
        nF   = size(F,1);
        
        nFacePerObject = (nF-nI)/nObjAdd;
        if ceil(nFacePerObject)-nFacePerObject > 0
            error('add:WrongInput','All patch objects have to be the same type!');
        end
        
        fIdxAdd = repmat([hAdd(:).number],[nFacePerObject 1]);
        fIdx = [fIdx; fIdxAdd(:)];
        setappdata(hPatch,'facenumber',fIdx);
        
        
    elseif strcmp(get(hPatch,'FaceColor'),'flat')
        fIdx = getappdata(hEPatch,'vertexnumber');
        nI = size(fIdx,1);
        V  = get(hEPatch,'Vertices');
        nV = size(V,1);
        
        nVertexPerObject = (nV-nI)/nObjAdd;
        if ceil(nVertexPerObject)-nVertexPerObject > 0
            error('add:WrongInput','All patch objects have to be the same type!');
        end
        
        fIdxAdd = repmat([hAdd(:).number],[nVertexPerObject 1]);
        fIdx = [fIdx; fIdxAdd(:)];
        setappdata(hPatch,'vertexnumber',fIdx);
        
    end
    
    if ~isempty(hTransform)
        set(hPatch,'Parent',hTransform);
    end
    set(hPatch,'Clipping','Off');
else
    % add all new objects to the hgtransform
    if ~isempty(hTransform)
        set(hNew,'Parent',hTransform);
    end
    set(hNew,'Clipping','Off');
end

% put lines to the bottom of the stack, otherwise shadows are
% wrong! MATLAB bug
if strcmp(get(hNew(1),'Tag'),'line')
    uistack(hNew,'bottom');
end

% comb together the handles of the old and new graphical objects.
sObject = [sObject hAdd(:)'];

% Saves the object handles into the figure UserData property.
setappdata(hFigure,'objects',sObject);
material(hAxis,'shiny');

% add callback function for showing the tooltips
%set([hNew hPatch],'ButtonDownFcn',@(obj,hit)swplot.tooltipcallback(obj,hit,hFigure,hTransform));
if showtooltip
    swplot.tooltip(swplot.tooltip,hFigure);
end

end