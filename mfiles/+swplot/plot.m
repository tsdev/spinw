function varargout = plot(varargin)
% plots objects to swplot figure
%
% SWPLOT.PLOT('Option1',Value1,...)
%
% hFigure = SWPLOT.PLOT(...)
%
%
% Options:
%
% type      Type of object to plot in a string. Possible options are:
%               'arrow'         position specifies start and end points
%               'ellipsoid'     position specifies center
%               'cylinder'      position specifies start and end points
%               'circle'        position specifies center and normal vector
%               'line'          position specifies start and end points
%               'text'          position specifies the center of the text
% position  Position of the object/objects in a matrix with dimensions of
%           [3 nObject 2]/[3 nObject] depending on the type of object.
% name      String, the name of the object. It can be used for finding the
%           object handles after plotting.
% text      Text to appear in the tooltip of the swplot figure after
%           clicking on the object. Can be a string that will be the same
%           for all objects, or a cell of strings for different text per
%           object. Default value is taken from the label option.
% label     Text to appear in the legend in a string for the same text of
%           all objects or strings in a cell for multiple objects with
%           dimension [1 nObject]. Default value is taken from the name
%           string.
% legend    Type of legend to show the object:
%               0       Do not show in legend.
%               1       Colored box in legend.
%               2       Dashed box in legend.
%               3       Colored sphere in legend.
% color     Color of objects, either a single color or as many colors as
%           many objects are given in a matrix with dimensions of [3 1]/[3
%           nObject]. Values are RGB triples with values between [0 255].
%           Can be also string or cell of strings with the name of the
%           colors, for details see swplot.color. Default is red.
% unit      String determining the coordinate system, either 'lu' for
%           lattice units where the lattice is defined by the stored
%           basis,or 'xyz' for the original matlab units. Default is 'lu'.
% figure    Handle of the swplot figure. Default is the selected figure.
% R         Radius value of cylinder, sphere (if no 'T' is given) and
%           arrow, default is 0.06.
% alpha     Head angle for arrow in degree units, default is 15 degree.
% lHead     Length of the arrow head, default value is 0.5.
% T         Transformation matrix that transforms a unit sphere to the
%           ellipse via: R' = T(:,:,i)*R
%           Dimensions are [3 3 nObject].
% nMesh     Resolution of the ellipse surface mesh. Integer number that is 
%           used to generate an icosahedron mesh with #mesh number of
%           additional triangulation, default value is stored in
%           swpref.getpref('nmesh')
% nPatch    Number of points on the curve for arrow and cylinder, default 
%           value is stored in swpref.getpref('npatch').
% hg        Whether to use hgtransform (nice rotation with the mouse) or 
%           default Matlab rotation of 3D objects. Default is true.
%
% See also SWPLOT.COLOR.
%

P0 = swpref.getpref('npatch',[]);
M0 = swpref.getpref('nmesh',[]);

inpForm.fname  = {'type' 'name' 'text' 'position' 'label' 'legend' 'color' 'unit' 'figure'};
inpForm.defval = {[]     []     []     []         []      []       []      'lu'   []      };
inpForm.size   = {[1 -8] [1 -1] [1 -2] [3 -3 -4]  [1 -5]  [1 1]    [-9 -6] [1 -7] [1 1]   };
inpForm.soft   = {false  true   true   false      true    true     true    false  true    };

inpForm.fname  = [inpForm.fname  {'R'     'alpha' 'lHead' 'nMesh' 'nPatch' 'T'       'hg'   }];
inpForm.defval = [inpForm.defval {0.06    15      0.5     M0      P0       []        'hg'   }];
inpForm.size   = [inpForm.size   {[1 -11] [1 1]   [1 1]   [1 1]   [1 1]    [3 3 -10] [1 -12]}];
inpForm.soft   = [inpForm.soft   {false   false   false   false   false    true      false}];

param = sw_readparam(inpForm, varargin{:});

type = lower(param.type);

if param.hg
    hgmode = 'hg';
else
    hgmode = 'nohg';
end

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
legend0 = [1 3 1 1 0 0];
type0   = {'arrow' 'ellipsoid' 'cylinder' 'circle' 'line'  'text' };
% default color per object type :D
col0    = {'red'   'blue'      'orange'   'gray'   'black' 'black'};
% create dictionary to convert string to number
K       = containers.Map(type0,1:6);
typeNum = K(type);

if isempty(param.legend)
    legend = legend0(typeNum);
else
    legend = param.legend;
end

if isempty(param.figure)
    hFigure = swplot.activefigure(hgmode);
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

% check size of position matrix
pos0 = [2 1 2 2 2 1];
pos  = param.position;
if size(pos,3) ~= pos0(typeNum)
    error('plot:WrongInput','The given position matrix has wrong dimensions!');
end

% number of objects to plot
nObject = size(pos,2);

% chekc unit selector string
switch lower(param.unit)
    case 'lu'
        BV = getappdata(hFigure,'base');
        % multiply the coordinates with the basis vectors
        xyz = permute(sum(bsxfun(@times,permute(pos,[4 1 2 3]),BV),2),[1 3 4 2]);
    case 'xyz'
        xyz = pos;
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

if nLabel ~= 1 && nLabel ~= nObject
    error('plot:WrongInput','Number of given labels does not agree with the number of object positions!')
end
if nCol ~= 1 && nCol ~= nObject
    error('plot:WrongInput','Number of given colors does not agree with the number of object positions!')
end


sObject = struct('handle',cell(1,nObject));

switch type
    case 'arrow'
        handle = swplot.arrow(hAxis,xyz(:,:,1),xyz(:,:,2),param.R,param.alpha,param.lHead,param.nPatch);
        
    case 'ellipsoid'
        if isempty(param.T)
            % use sphere drawing mode generating spheres from option 'R'
            param.T = bsxfun(@times,eye(3),permute(param.R,[1 3 2]));
            if numel(param.T) == 9
                param.T = repmat(param.T,[1 1 nObject]);
            end
        end
        
        handle = swplot.ellipsoid(hAxis,xyz,param.T,param.nMesh);
    case 'cylinder'
        % closed cylinder
        handle = swplot.cylinder(hAxis,xyz(:,:,1),xyz(:,:,2),param.R,param.nPatch,true);
    case 'circle'
        handle = gobjects(1,nObject);
        for ii = 1:nObject
            handle(ii) = swplot.circle(hAxis,xyz(:,ii,1),xyz(:,ii,2),param.R,param.nPatch);
        end
        % remove normal vectors (not stored in appdata)
        pos = pos(:,1,:);
    case 'line'
        handle = gobjects(1,nObject);
        for ii = 1:nObject
            handle(ii) = swplot.line(hAxis,xyz(:,ii,1),xyz(:,ii,2));
        end
    case 'text'
        textStr = param.text;
        if ~iscell(textStr)
            textStr = repmat({textStr},[1 nObject]);
        end
        
        handle = swplot.text(hAxis,xyz,textStr);
end

handle = num2cell(handle);
[sObject(:).handle] = handle{:};

% change color of objects
if sObject(1).handle == getappdata(hFigure,'facepatch')
    % there is only a single unique handle
    hPatch = sObject(1).handle;
    
    if nCol == 1
        patchCData = repmat(color,[1 nObject]);
    else
        patchCData = color;
    end
    
    % set colors per face
    fIdx = getappdata(hPatch,'facenumber');
    nI   = size(fIdx,1);
    C    = get(hPatch,'FaceVertexCData');
    nC   = size(C,1);
    
    nNewFace = nC-nI;
    nFacePerObject = nNewFace/nObject;
    if ceil(nFacePerObject)-nFacePerObject > 0
        error('plot:WrongInput','All patch objects have to be the same type!');
    end
    
    patchCData = reshape(permute(repmat(patchCData,[1 1 nFacePerObject]),[3 2 1]),[],3);
    C(end+(((-nNewFace+1):0)),:) = patchCData;
    set(hPatch,'FaceVertexCData',C);

else
    % set color per handle
    % change color of object using the right property prop0
    if strcmp(get(sObject(1).handle,'type'),'patch')
        prop0 = 'FaceColor';
    else
        prop0 = 'Color';
    end
    
    if nCol == 1
        % same color for every object
        set([sObject(:).handle],prop0,color);
    else
        % different color for each object
        for ii = 1:nObject
            set([sObject(ii).handle],prop0,color(:,ii));
        end
    end
end

% save text
if ~iscell(param.text)
    param.text = repmat({param.text},[1 nObject]);
end
[sObject(:).text] = param.text{:};

% save position
posC = mat2cell(permute(pos,[1 3 2]),3,size(pos,3),ones(1,nObject));
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

% add objects to the figure
swplot.add(sObject,hFigure);

% take care of the legend
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

if nargout > 0
    varargout{1} = hFigure;
end

% final corrections
swplot.zoom

end