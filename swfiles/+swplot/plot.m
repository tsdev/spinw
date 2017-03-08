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
%               'polyhedron'    position specifies the vertices of the
%                               convex polyhedron or polygon
%               'circle'        position specifies center and normal vector
%               'line'          position specifies start and end points (or
%                               any number of points per curve)
%               'text'          position specifies the center of the text
%               'orbital'       position specifies the center of the
%                               orbital
% position  Position of the object/objects in a matrix with dimensions of
%           [3 nObject 2]/[3 nObject]/[3 nObject nPoint] depending on the
%           type of object.
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
%           colors, for details see swplot.color. Default is red. For
%           orbitals two colors per orbital can be defined, for the lobes
%           corresponding to the +/- values of the wave function. The
%           second set of colors can be added along the third dimension,
%           e.g. [3 1 2]/[3 nObject 2] for RGB values and string-cell with
%           dimensions of [1 1 2]/[1 nObject 2]. Default value is red/blue.
% alpha     Transparency of objects (1 non-transparent, 0 transparent)
%           defined as a single number for unitform transparency or as a
%           row vector with nObject element to set transparency per object.
%           Default value is 1.
% unit      String that determines the coordinate system:
%               'lu'    Lattice units are used where the lattice is defined
%                       by the stored basis (default).
%               'xyz'   Use the original matlab units.
% figure    Handle of the swplot figure. Default is the selected figure.
% R         Radius value of cylinder, sphere (if no 'T' is given) and
%           arrow, default is 0.06.
% ang       Angle for arrow head in degree units, default is 15 degree.
% lHead     Length of the arrow head, default value is 0.5.
% T         Transformation matrix that transforms objects:
%               'ellipsoid'     a unit sphere is transformed to an ellipsoid:
%                                   R' = T(:,:,i)*R
%               'orbital'       a Hydrogen orbital is rotated to arbitrary
%                               orientation:
%                                   O' = T(:,:,i)*O*T(:,:,i)'
%           Dimensions are [3 3 nObject].
% lineStyle Line style, default value is '-' for continuous lines. It can
%           be also a vector with as many elements as many line segments.
%           In this case the numbers are equivalent to the following style
%           format string:
%               1   '-'
%               2   '--'
%               3   '-.'
%               4   ':'
%               5   'none'
% lineWidth Line width, default value is 0.5, can be a vector with nObject
%           columns for different width per line segment.
% scale     Scale for the orbitals, the default scale=1 corresponds to the
%           true size of the Hydrogen orbital in Angstrom.
% qLabel    Label of the orbital in a string, for details check the same
%           option in swplot.orbital.
% nMesh     Resolution of the ellipse surface mesh. Integer number that is
%           used to generate an icosahedron mesh with #mesh number of
%           additional triangulation, default value is stored in
%           swpref.getpref('nmesh')
% nPatch    Number of points on the curve for arrow and cylinder, default
%           value is stored in swpref.getpref('npatch').
% tooltip   If true, the tooltip will be switched on at the end of the
%           plot. Default is true.
% replace   If true, all object with the same name as the new plot will be
%           deleted before plotting. Default is false.
% data      Arbitrary data per object that will be stored in the swplot
%           figure and can be retrieved. It is stored in a cell with
%           nObject number of elements.
% translate If true, all plot objects will be translated to the figure
%           center. Default is true.
% zoom      If true, figure will be automatically zoomed to the ideal size.
%           Default is true.
%
% See also SWPLOT.COLOR.
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

inpForm.fname  = [inpForm.fname  {'scale' 'qLabel'}];
inpForm.defval = [inpForm.defval {1       ''      }];
inpForm.size   = [inpForm.size   {[1 1]   [1 -18] }];
inpForm.soft   = [inpForm.soft   {false   true    }];

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
% 7 'orbital'
legend0 = [1 3 1 1 0 0 1 0];
type0   = {'arrow' 'ellipsoid' 'cylinder' 'circle' 'line'  'text' 'polyhedron' 'orbital'};
% default color per object type :D
col0    = {'red'   'blue'      'orange'   'gray'   'black' 'black' 'turquoise' cat(3,{'red'},{'blue'})};
% create dictionary to convert string to number
K       = containers.Map(type0,1:numel(col0));
typeNum = K(type);

if strcmp(type,'orbital') && isempty(param.qLabel)
    error('plot:OrbitalMissing','Missing qLabel option!')
end

if isempty(param.legend)
    legend = legend0(typeNum);
else
    legend = param.legend;
end

if isempty(param.name)
    param.name = type;
end

if isempty(param.label)
    param.label = param.name;
end

% label for legend
if ~iscell(param.label)
    label = {param.label};
else
    label = param.label;
end
nLabel = numel(label);

% check size of position matrix
pos0 = [2 1 2 2 0 1 0 1];
pos  = param.position;
if pos0(typeNum)~=0 && size(pos,3) ~= pos0(typeNum)
    error('plot:WrongInput','The given position matrix has wrong dimensions!');
end

% number of objects to plot
nObject = size(pos,2);

% check colors
if iscellstr(param.color) && size(param.color,1)~=1
    error('color:WrongInput','Color option has wrong dimensions!');
end
if size(param.color,3)>1 && ~strcmp(type,'orbital')
    error('color:WrongInput','Color option has wrong dimensions!')
end

if isempty(param.color)
    param.color = col0{typeNum};
end

% convert colors to Matlab color values
if ischar(param.color)
    sColor = [1 1];
else
    sColor = size(param.color);
end

if isnumeric(param.color)
    color  = swplot.color(reshape(param.color,3,[]))/255;
else
    color  = swplot.color(param.color)/255;
end
sColor(1) = 3;
color  = reshape(color,sColor);

% number of colors
nCol  = size(color,2);

% correct color dimensions for orbital type plot
if strcmp(type,'orbital')
    if nCol==2 && nCol~=nObject && size(color,3)==1
        color = permute(color,[1 3 2]);
        nCol  = 1;
    end
end

% number of transparency values
nAlp  = numel(param.alpha);

if nLabel ~= 1 && nLabel ~= nObject
    error('plot:WrongInput','Number of given labels does not agree with the number of object positions!')
end
if ~strcmp(type,'orbital') && nCol ~= 1 && nCol ~= nObject
    error('plot:WrongInput','Number of given colors does not agree with the number of object positions!')
end
if strcmp(type,'orbital') && nCol ~= 1 && nCol ~= 2 && nCol ~= nObject && nCol ~= 2*nObject
    error('plot:WrongInput','Number of given colors does not agree with the number of object positions!')
end

if isempty(param.figure)
    hFigure = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

% axis handle
hAxis = getappdata(hFigure,'axis');

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

% remove objects with the same name if necessary
if param.replace
    % find objects to be deleted
    sObj = swplot.findobj(hFigure,'name',param.name);
    % delete them!
    swplot.delete([sObj(:).number]);
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
    case 'orbital'
        if isempty(param.T)
            param.T = eye(3);
        end
        
        if numel(param.T) == 9
            param.T = repmat(param.T,[1 1 nObject]);
        end
        
        handle = swplot.orbital(hAxis,param.qLabel,xyz,param.T,param.scale,param.nPatch);
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

        if strcmp(type,'orbital')
            % allow 2 colors per orbital: 
            %   color1: (:,1:nObject)
            %   color2: (:,(nObject+1):end)
            if size(color,3) == 2
                % different colors for the +/- lobes
                patchCData = [repmat(color(:,:,1),[1 nObject/nCol]) repmat(color(:,:,2),[1 nObject/nCol])];
            else
                % same color for the +/- lobes
                patchCData = repmat(color,[1 2*nObject/nCol]);
            end
        else
            patchCData = repmat(color,[1 nObject/nCol]);
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
        
        if strcmp(type,'orbital')
            % assign the two colors
            patchCData1 = reshape(permute(repmat(patchCData(:,1:(end/2)),[1 1 nFacePerObject]),[3 2 1]),[],3);
            patchCData2 = reshape(permute(repmat(patchCData(:,(end/2+1):end),[1 1 nFacePerObject]),[3 2 1]),[],3);
            % mix in CData2 where C is blue C(:,3)
            blueIdx = logical(C(end+(((-nNewFace+1):0)),3));
            patchCData1(blueIdx,:) = patchCData2(blueIdx,:);
            C(end+(((-nNewFace+1):0)),:) = patchCData1;
        else
            patchCData = reshape(permute(repmat(patchCData,[1 1 nFacePerObject]),[3 2 1]),[],3);
            C(end+(((-nNewFace+1):0)),:) = patchCData;
        end
        
        patchAlphaData = reshape(repmat(patchAlphaData,[nFacePerObject 1]),[],1);
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
        % same color for every object except orbitals
        if strcmp(type,'orbital')
            nFace = size(get(sObject(1).handle,'Faces'),1);
            C1 = repmat(color(:,1,1),[1 nFace 1])';
            C2 = repmat(color(:,1,2),[1 nFace 1])';
            C0 = get(sObject(1).handle,'FaceVertexCData');
            blueIdx = logical(C0(:,3));
            C1(blueIdx,:) = C2(blueIdx,:);
            set([sObject(:).handle],'FaceVertexCData',C1);
        else
            set([sObject(:).handle],propC,color);
        end
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