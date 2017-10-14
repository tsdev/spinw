function varargout = plot(varargin)
% plots objects to swplot figure
% 
% ### Syntax
% 
% `swplot.plot(Name,Value)`
%
% `hFigure = swplot.plot(Name,Value)`
% 
% ### Description
% 
% `swplot.plot(Name,Value)` plots objects to the swplot figure and adds the
% objects to the [matlab.hgtransform] object. This command enables the
% plotting of multiple objects simultaneously while enabling fine control
% of color, legend test, tooltip text etc. This commands is used by the
% [spinw.plot] high level plot command.
% 
% ### Name-Value Pair Arguments
% 
% `'type'`
% : Type of object to plot in a string. Possible options are:
%   * `'arrow'`         position specifies start and end points,
%   * `'ellipsoid'`     position specifies center,
%   * `'cylinder'`      position specifies start and end points,
%   * `'polyhedron'`    position specifies the vertices of the
%                       convex polyhedron or polygon,
%   * `'circle'`        position specifies center and normal vector,
%   * `'line'`          position specifies start and end points (or
%                       any number of points per curve),
%   * `'text'`          position specifies the center of the text.
% 
% `'position'`
% : Position of the object/objects in a matrix with dimensions of
%   $[3\times n_{obj}\times 2]$ or $[3\times n_{obj}]$ or $[3\times
%   n_{obj}\times n_{point}]$ depending on the type of object. The unit of
%   the positions is determined by the `unit` parameter.
% 
% `'name'`
% : String, the name of the object. It can be used for grouping the
%   object handles to enable easier search, see [swplot.findobj] for
%   details.
% 
% `'text'`
% : Text to appear in the tooltip of the swplot figure after
%   clicking on the object. Can be a string that will be the same
%   for all objects, or a cell of strings for different text per
%   object. Default value is taken from the label option.
% 
% `'label'`
% : Text to appear in the legend in a string for the same text of
%   all objects or strings in a cell with $n_{obj}$ number of elements for
%   multiple objects. Default value is taken from the `name` parameter.
% 
% `'legend'`
% : Type of legend to show the object:
%   * `0`       do not show in legend,
%   * `1`       colored box in legend,
%   * `2`       dashed box in legend,
%   * `3`       colored sphere in legend.
% 
% `'color'`
% : Color of objects, either a single color or as many colors as
%   many objects are given in a matrix with dimensions of $[3\times 1]$ or
%   $[3\times n_{obj}]$. Colors are RGB triplets with values between 0 and
%   255. Can be also string or cell of strings with the name of the colors,
%   for possible color names see [swplot.color]. Default value is `'red'`.
% 
% `'alpha'`
% : Transparency of objects (1: non-transparent, 0: transparent)
%   defined as a single number for uniform transparency or as a
%   row vector with $n_{obj}$ number of elements to set transparency per object.
%   Default value is 1.
% 
% `'unit'`
% : String that determines the coordinate system where position vectors are
%   defined:
%   * `'lu'`    Lattice units are used where the lattice is defined
%               by the stored basis (default).
%   * `'xyz'`   Use the original Matlab units.
% 
% `'figure'`
% : Handle of the swplot figure, default is the active figure.
% 
% `'R'`
% : Radius value of cylinder, sphere (if no `'T'` parameter is given) and
%   arrow, default value is 0.06.
% 
% `'ang'`
% : Angle for arrow head in degree, default value is 15\\deg.
% 
% `'lHead'`
% : Length of the arrow head, default value is 0.5.
% 
% `'T'`
% : Transformation matrix that transforms a unit sphere to the
%   ellipse via: `R' = T(:,:,i)*R`, stored in a matrix with
%   dimensions of $[3\times 3\times n_{obj}]$.
% 
% `'lineStyle'`
% : Line style, default value is `'-'` for continuous lines. It can
%   be also a vector with as many elements as many line segments.
%   In this case the numbers are equivalent to the following style
%   format string:
%   * `1`   `'-'`,
%   * `2`   `'--'`,
%   * `3`   `'-.'`,
%   * `4`   `':'`,
%   * `5`   `'none'`.
% 
% `'lineWidth'`
% : Line width, default value is 0.5, can be a vector with $n_{obj}$
%   columns for different width per line segment.
% 
% `'fontSize'`
% : Font size of text in pt when `type` parameter is set to `'text'`.
%   Default value is stored in `swpref.getpref('fontsize')`.
% 
% `'nMesh'`
% : Resolution of the ellipse surface mesh. Integer number that is
%   used to generate an icosahedron mesh with `nMesh` number of
%   additional subdivision of triangular surfaces. Default value is stored in
%   `swpref.getpref('nmesh')`.
% 
% `'nPatch'`
% : Number of points on the curve for arrow and cylinder, default
%   value is stored in `swpref.getpref('npatch')`.
% 
% `'tooltip'`
% : If `true`, the tooltip will be switched on after the
%   plot. Default is `true`.
% 
% `'replace'`
% : If `true`, all objects with the same name as the new plot will be
%   deleted before plotting. Default is `false`.
% 
% `'data'`
% : User supplied data per object that will be stored in the swplot
%   figure and can be retrieved using [swplot.getdata]. It is stored in a
%   cell with $n_{obj}$ number of elements.
% 
% `'translate'`
% : If `true`, the average center of the plot objects will be translated to
%   the figure center. Default is `true`.
% 
% `'zoom'`
% : If `true`, the swplot figure will be zoomed to make the plot objects
%   cover the full figure. Default is `true`.
% 
% ### See Also
% 
% [swplot.color] \| [swplot.add]
%

P0 = swpref.getpref('npatch',[]);
M0 = swpref.getpref('nmesh',[]);
fontSize0 = swpref.getpref('fontsize',[]);

inpForm.fname  = {'type' 'name' 'text' 'position' 'label' 'legend' 'color' 'unit' 'figure' 'lineStyle'};
inpForm.defval = {[]     []     ''     []         []      []       []      'lu'   []       '-'        };
inpForm.size   = {[1 -8] [1 -1] [1 -2] [3 -3 -4]  [1 -5]  [1 1]    [-9 -6] [1 -7] [1 1]   [1 -13]     };
inpForm.soft   = {false  true   true   false      true    true     true    false  true    false       };

inpForm.fname  = [inpForm.fname  {'R'     'ang'   'lHead' 'nMesh' 'nPatch' 'T'       'hg'    'tooltip'}];
inpForm.defval = [inpForm.defval {0.06    15      0.5     M0      P0       []        'hg'    true     }];
inpForm.size   = [inpForm.size   {[1 -11] [1 1]   [1 1]   [1 1]   [1 1]    [3 3 -10] [1 -12] [1 1]    }];
inpForm.soft   = [inpForm.soft   {false   false   false   false   false    true      false   false    }];

inpForm.fname  = [inpForm.fname  {'data'    'replace' 'lineWidth' 'fontSize' 'translate' 'zoom' 'alpha'}];
inpForm.defval = [inpForm.defval {{}        false     0.5         fontSize0  true        true   1      }];
inpForm.size   = [inpForm.size   {[-14 -15] [1 1]     [1 -16]     [1 1]      [1 1]       [1 1]  [1 -17]}];
inpForm.soft   = [inpForm.soft   {true      false     false       false      false       false  false  }];

param = sw_readparam(inpForm, varargin{:});

type = lower(param.type);

% convert type string to index
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
type0   = {'arrow' 'ellipsoid' 'cylinder' 'circle' 'line'  'text' 'polyhedron'};
% default color per object type :D
col0    = {'red'   'blue'      'orange'   'gray'   'black' 'black' 'turquoise'};
% create dictionary to convert string to number
K       = containers.Map(type0,1:numel(col0));
typeNum = K(type);

if isempty(param.legend)
    legend = legend0(typeNum);
else
    legend = param.legend;
end

if isempty(param.figure)
    hFigure = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

% axis handle
hAxis = getappdata(hFigure,'axis');

if isempty(param.name)
    param.name = type;
end

if isempty(param.label)
    param.label = param.name;
end

% if isempty(param.text)
%     param.text = param.label;
% end

% label for legend
if ~iscell(param.label)
    label = {param.label};
else
    label = param.label;
end
nLabel = numel(label);

% remove objects with the same name if necessary
if param.replace
    % find objects to be deleted
    sObj = swplot.findobj(hFigure,'name',param.name);
    % delete them!
    swplot.delete([sObj(:).number]);
end

% check size of position matrix
pos0 = [2 1 2 2 0 1 0];
pos  = param.position;
if pos0(typeNum)~=0 && size(pos,3) ~= pos0(typeNum)
    error('plot:WrongInput','The given position matrix has wrong dimensions!');
end

% number of objects to plot
nObject = size(pos,2);

% basis vectors
BV = getappdata(hFigure,'base');

% check unit selector string
% pos will store coordinates in lu units
switch lower(param.unit)
    case 'lu'
        % multiply the coordinates with the basis vectors
        xyz = permute(sum(bsxfun(@times,permute(pos,[4 1 2 3]),BV),2),[1 3 4 2]);
    case 'xyz'
        
        xyz = pos;
        pos = permute(sum(bsxfun(@times,permute(pos,[4 1 2 3]),inv(BV)),2),[1 3 4 2]);
    otherwise
        error('plot:WrongInput','The selected coordinate system unit option does not exists!');
end

% check colors
if iscellstr(param.color) && size(param.color,1)~=1
    error('color:WrongInput','Color option has wrong dimensions!');
end

if isempty(param.color)
    param.color = col0{typeNum};
end

% convert colors to Matlab color values
color = swplot.color(param.color)/255;
% number of colors
nCol  = size(color,2);
% number of transparency values
nAlp  = numel(param.alpha);

if nLabel ~= 1 && nLabel ~= nObject
    error('plot:WrongInput','Number of given labels does not agree with the number of object positions!')
end
if nCol ~= 1 && nCol ~= nObject
    error('plot:WrongInput','Number of given colors does not agree with the number of object positions!')
end


sObject = struct('handle',cell(1,nObject));

switch type
    case 'arrow'
        handle = swplot.arrow(hAxis,xyz(:,:,1),xyz(:,:,2),param.R,param.ang,param.lHead,param.nPatch);
        
    case 'ellipsoid'
        if isempty(param.T)
            % use sphere drawing mode generating spheres from option 'R'
            param.T = bsxfun(@times,eye(3),permute(param.R,[1 3 2]));
            if numel(param.T) == 9
                param.T = repmat(param.T,[1 1 nObject]);
            end
        end
        
        handle = swplot.ellipsoid(hAxis,xyz,param.T,param.nMesh);
        pos(:,:,2) = nan;
    case 'polyhedron'
        handle = swplot.polyhedron(hAxis,xyz);
        % get the extremum positions
        pos = cat(3,min(pos,[],3),max(pos,[],3));

    case 'cylinder'
        % closed cylinder
        handle = swplot.cylinder(hAxis,xyz(:,:,1),xyz(:,:,2),param.R,param.nPatch,true);
    case 'circle'
        handle = swplot.circle(hAxis,xyz(:,:,1),xyz(:,:,2),param.R,param.nPatch);
        % remove normal vectors (use nans)
        pos(:,:,2) = nan;
    case 'line'
        handle = swplot.line(hAxis,xyz,[],param.lineStyle,param.lineWidth,true);
        % get the extremum positions
        pos = cat(3,min(pos,[],3),max(pos,[],3));
    case 'text'
        textStr = param.text;
        if ~iscell(textStr)
            textStr = repmat({textStr},[1 nObject]);
        end
        pos(:,:,2) = nan;
        handle = swplot.text(hAxis,xyz,textStr,param.fontSize);
end

handle = num2cell(handle);
[sObject(:).handle] = handle{:};

% change color and transparency of objects
if numel(handle)>1 && handle{1}==handle{2}
    % all handle point to the same object, there is only a single unique
    % handle
    hPatch = sObject(1).handle;
    
    if strcmp(get(hPatch,'FaceColor'),'flat')
        % it is a faces object
        
        if nCol == 1
            patchCData = repmat(color,[1 nObject]);
        else
            patchCData = color;
        end
        if nAlp == 1
            patchAlphaData = repmat(param.alpha,[1 nObject]);
        else
            patchAlphaData = param.alpha;
        end
        
        % set colors per face
        fIdx = getappdata(hPatch,'facenumber');
        nI   = size(fIdx,1);
        C    = get(hPatch,'FaceVertexCData');
        A    = get(hPatch,'FaceVertexAlphaData');
        nC   = size(C,1);
        
        nNewFace = nC-nI;
        nFacePerObject = nNewFace/nObject;
        if ceil(nFacePerObject)-nFacePerObject > 0
            error('plot:WrongInput','All patch objects have to be the same type!');
        end
        
        patchCData = reshape(permute(repmat(patchCData,[1 1 nFacePerObject]),[3 2 1]),[],3);
        patchAlphaData = reshape(repmat(patchAlphaData,[nFacePerObject 1]),[],1);
        C(end+(((-nNewFace+1):0)),:) = patchCData;
        A(end+(((-nNewFace+1):0)),:) = patchAlphaData;
        set(hPatch,'FaceVertexCData',C,'FaceVertexAlphaData',A);
        
    elseif strcmp(get(hPatch,'EdgeColor'),'flat')
        % object has only edges
        
        if nCol == 1
            patchCData = repmat(color,[1 nObject]);
        else
            patchCData = color;
        end
        if nAlp == 1
            patchAlphaData = repmat(param.alpha,[1 nObject]);
        else
            patchAlphaData = param.alpha;
        end
        
        % set colors per edge
        fIdx = getappdata(hPatch,'vertexnumber');
        nI   = size(fIdx,1);
        C    = get(hPatch,'FaceVertexCData');
        A    = get(hPatch,'FaceVertexAlphaData');
        nC   = size(C,1);
        
        nNewEdge = nC-nI;
        nEdgePerObject = nNewEdge/nObject;
        if ceil(nEdgePerObject)-nEdgePerObject > 0
            error('plot:WrongInput','All patch objects have to be the same type!');
        end
        
        patchCData = reshape(permute(repmat(patchCData,[1 1 nEdgePerObject]),[3 2 1]),[],3);
        patchAlphaData = reshape(repmat(patchAlphaData,[nEdgePerObject 1]),[],1);
        C(end+(((-nNewEdge+1):0)),:) = patchCData;
        A(end+(((-nNewEdge+1):0)),:) = patchAlphaData;
        set(hPatch,'FaceVertexCData',C,'FaceVertexAlphaData',A);
    else
        error('plot:WrongPatchObject','Patch object is wrong!');
    end
else
    % set color and transparency per handle
    % change color of object using the right property prop0
    if strcmp(get(sObject(1).handle,'type'),'patch')
        if strcmp(get(sObject(1).handle,'FaceColor'),'none')
            propC = 'EdgeColor';
            propA = 'EdgeAlpha';
        else
            propC = 'FaceColor';
            propA = 'FaceAlpha';
        end
    else
        propC = 'Color';
        propA = '';
    end
    
    if nCol == 1
        % same color for every object
        set([sObject(:).handle],propC,color);
    else
        % different color for each object
        for ii = 1:nObject
            set([sObject(ii).handle],propC,color(:,ii));
        end
    end
    if ~isempty(propA)
        if nAlp == 1
            % same transparency for every object
            set([sObject(:).handle],propA,param.alpha);
        else
            % different transparency for each object
            for ii = 1:nObject
                set([sObject(ii).handle],propA,param.alpha(ii));
            end
        end
    end
    
end

% save text
if ~iscell(param.text)
    param.text = repmat({param.text},[1 nObject]);
end
[sObject(:).text] = param.text{:};

% save position
if nObject == 1
    posC = mat2cell(permute(pos,[1 3 2]),3,size(pos,3));
else
    posC = mat2cell(permute(pos,[1 3 2]),3,size(pos,3),ones(1,nObject));
end
[sObject(:).position] = posC{:};

% save label
if nLabel == 1
    labelC = repmat(label,[1 nObject]);
else
    labelC = label;
end
[sObject(:).label] = labelC{:};

% save legend
legendC = repmat({param.legend},[1 nObject]);
[sObject(:).legend] = legendC{:};

% save type
typeC = repmat({type},[1 nObject]);
[sObject(:).type] = typeC{:};

% save name
nameC = repmat({param.name},[1 nObject]);
[sObject(:).name] = nameC{:};

% save data
if ~isempty(param.data)
    [sObject(:).data] = param.data{:};
end

% add objects to the figure
swplot.add(sObject,hFigure,param.tooltip);

% take care of the legend
if param.legend
    lDat = getappdata(hFigure,'legend');
    if nLabel > 1 || nCol > 1
        % different legend per object
        if nCol == 1
            color = repmat(color,[1 nObject]);
        end
        if nLabel == 1
            label = repmat(label,[1 nObject]);
        end
        
        legend = repmat(legend,[1 nObject]);
    end
    
    % append text
    if ~isempty(lDat.text)
        lDat.text = [lDat.text(:)' label];
    else
        lDat.text = label;
    end
    % append color
    lDat.color = [lDat.color color];
    % append type
    lDat.type = [lDat.type legend];
    
    setappdata(hFigure,'legend',lDat);
end

if nargout > 0
    varargout{1} = hFigure;
end

% final corrections to make figure nice
if param.zoom
    swplot.zoom('auto',hFigure);
end

if param.translate
    swplot.translate('auto',hFigure);
end

end