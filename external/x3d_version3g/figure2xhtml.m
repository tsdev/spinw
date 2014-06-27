function figure2xhtml(varargin)
% This function FIGURE2XHTML converts the 3D objects of a Matlab figure
% to an in XHTML embedded X3D file. In a modern browser, you can now view
% your figure interactively because of the great X3DOM Library
% (Instant 3D the HTML way! http://www.x3dom.org/).
%
% Currently the function supports: Axes, Patch, Line, Surface, Text, Image
% and Light objects.
%
% Browsers supported:
% - Google Chrome 9.x and above
% - Firefox 4.x and above
% - Webkit nightly builds
% - InstantReality-plugin or Flash11-plugin needed for Internet explorer
%
%
%  figure2xhtml(filename,handle,options)
%
%   Note : All input arguments are optional, and sorting doesn't matter
%          thus figure2xhtml(); is valid but also
%          figurexhtml(options,handle);
%
% inputs,
%    filename : Name of the XHTML file (also an X3D file is created).
%           When empty or not known a File-Dialog is shown
%    handle : Figure handle or axis handle
%           When empty or not known the current axis (GCA) is used
%    options : A struct with options
%      options.output : Produce output files 'x3d', 'xhtml' or 'both' (default)
%      options.width : Width of X3D render object in pixels default 500 
%      options.height : Height of X3D render object in pixels default 500 
%      options.headlight : Enable Camera head light, boolean true/false
%                           (default true)
%      options.embedimages : Instead of using separate .png files embed 
%                            the images in xhtml (Not supported bij IE9)
%                           (default false)
%      options.title : Title of xhtml page, default 'Matlab X3D'
%      options.interactive : Make mesh/surface objects clickable in xhtml,
%                          boolean true/false (default false)
%                           
%
% Example, Click-able patch, with face color
%  load('functions\exampledata');
%  figure, hold on; axis equal;
%  patch(FV,'facecolor',[1 0 0],'facealpha',0.5);
%  FV.vertices(:,1)=FV.vertices(:,1)+80;
%  patch(FV,'FaceColor','interp','FaceVertexCData',rand(size(FV.vertices)),'edgecolor','none');
%  FV.vertices(:,1)=FV.vertices(:,1)+80;
%  patch(FV,'FaceColor','flat','FaceVertexCData',rand(size(FV.faces)),'edgecolor','none');
%  figure2xhtml('test/example1',struct('interactive',true))
%
%Example, Lights and surface
%   logo
%   if(exist('l1','var'))
%     set(l1,'Position',[40 100 20]*2);
%     set(l2,'Position',[.5 -1 .4]*30);
%   end
%   figure2xhtml('test/example2.xhtml')
%
%Example, Surf with colors
%   h=figure, hold on; sphere;
%   figure2xhtml('test/example3',h)
%
%Example, Mesh
%   figure, axis([-3 3 -3 3 -10 5])
%   [X,Y] = meshgrid(-3:.125:3);
%   Z = peaks(X,Y);
%   mesh(X,Y,Z);
%   figure2xhtml('test/example4')
%
%Example, Surf
%   figure,
%   [X,Y] = meshgrid(-3:.125:3);
%   Z = peaks(X,Y);
%   surf(X,Y,Z);
%   figure2xhtml;
%
%Example, Lines
%   figure, hold on;
%   plot3([1 10 30]*4,[1 20 30]*4,[1 10 30]*4,'r:');
%   plot3([1 10 30]'*4,[1 20 30]'*4,[10 30 1]'*4,'g*');
%   plot3([1 5 20]'*4,[1 5 20]'*4,[3 20 1]'*4,'-g*','MarkerEdgeColor',[1 0 0]);
%   figure2xhtml('test/example5')
%
%Example, Points
%   L = 40*membrane(1,25);
%   [X,Y]=ndgrid(1:size(L,1),1:size(L,2));
%   figure, plot3(X(:),Y(:),L(:),'r.');
%   figure2xhtml('test/example6')
%
%Example, Texture on Surface
%   figure, hold on;
%   load('functions\exampledata');
%   surface(peaks,flipud(X),'FaceColor','texturemap','EdgeColor','none', 'CDataMapping','direct')
%   colormap(map)
%   view(-35,45)
%   figure2xhtml('test/example7')
%
%Example, Axis
%   h=figure, hold on; sphere;
%   addpath([cd '/functions'])
%   axis2lines
%   figure2xhtml('test/example8',h)
%
% Function is written by D.Kroon University of Twente (July 2011)

% Check input arguments
haxis=[]; filename=[]; options=[];
if(nargin>3)
   error('figure2xhtml:inputs','too much input variables');
end
for ii = 1:nargin
    val = varargin{ii};
    if(ischar(val)),
        filename=val;
    elseif(isnumeric(val)),
        haxis=val;
    elseif(isstruct(val)), 
        options=val;
    elseif isa(gcf,'matlab.ui.Figure')
        haxis = val;
    else
        error('figure2xhtml:inputs','unknown input');
    end
end

% Process inputs
defaultoptions=struct('output','both','height',500,'width',500,'headlight',true,'title','Matlab X3D','interactive',false,'embedimages',false);
if(isempty(options)), options=defaultoptions;
else
    tags = fieldnames(defaultoptions);
    for ii=1:length(tags), if(~isfield(options,tags{ii})), options.(tags{ii})=defaultoptions.(tags{ii}); end, end
    if(length(tags)~=length(fieldnames(options))),
        warning('register_images:unknownoption','unknown options found');
    end
end

% Get axis handle of figure-handle or current handle
haxis=axis_handle(haxis);

if(isempty(filename))
    [filename, folder] = uiputfile({'*.xhtml','XHTML file (*.xhtml)';'*.x3d','X3D (*.x3d)'; '*.*', 'All Files (*.*)'}, 'Select an output filename');
    filename=[folder filename];
end

% Extension  of filename
[folder,filename,ext]=fileparts(filename);
if(~isempty(folder)), folder=[folder '/']; end

% Causes all graphics objects to be ready
drawnow('expose'); pause(0.1);

% add all needed function paths
add_function_paths;

% get Tag properties
Ch = get(haxis,'Children');
Tag0 = {};

for ii = 1:numel(Ch)
    switch(Ch(ii).Type)
        case 'patch'
            Tag0{end+1} = get(Ch(ii),'Tag'); %#ok<*AGROW>
        case 'line'
        case 'surface'
            Tag0{end+1} = get(Ch(ii),'Tag');
        case 'image'
            Tag0{end+1} = get(Ch(ii),'Tag');
        case 'light'
        case 'text'
    end
end

Tag0 = regexprep(Tag0,'_',': ');

% Create the xhtml file
if(strcmpi(options.output,'xhtml')||strcmpi(options.output,'both'))
    [data,loc_body]=XHTMLheader(options);
    data.tags.numobjects=0;
    data.tags.xhtml=true;
    data.tags.folder=folder;
    data.tags.filename=filename;
    data.tags.options=options;
    [data,loc_scene]=X3Dheader(data,loc_body,haxis);
    data=figurex3d(haxis,data,loc_scene);
    if(options.interactive)
        data=XMLaddNode('script',data,loc_body+1);
        data=XMLaddProperty('type','text/javascript',data);
        data=XMLaddString(interactive_js(data.tags.numobjects,Tag0),data);
    end
    writexhtmlfile(folder,filename,data);
end

% Create the x3d file
if(strcmpi(options.output,'x3d')||strcmpi(options.output,'both'))
    data=struct;
    data.tags.numobjects=0;
    data.tags.xhtml=false;
    data.tags.folder=folder;
    data.tags.filename=filename;
    data.tags.options=options;
    [data,loc_scene]=X3Dheader(data,[],haxis);
    data=figurex3d(haxis,data,loc_scene);
    writex3dfile(folder,filename,data);
end

function data=figurex3d(handle,data,loc_scene)
% Get all Childeren of the current axis, and create an X3D xml structure
% of these 3D objects
Ch=get(handle,'Children');
for i=1:length(Ch)
    Obj=get(Ch(i));
    type=Obj.Type;
    switch(type)
        case 'axes'
        case 'patch'
            data=addmesh(data,loc_scene,Obj);
        case 'line'
            data=addline(data,loc_scene,Obj);
        case 'surface'
            [F,V,Cface,Cedge,E,T]=surf2FV(Obj);
            Obj.Faces=F;
            Obj.E=E;
            Obj.Vertices=V;
            Obj.FaceVertexCData=Cface;
            Obj.EdgeVertexCData=Cedge;
            Obj.TextureVertices=T;
            data=addmesh(data,loc_scene,Obj);
        case 'image'
            C=Obj.CData;
            xd=[Obj.XData(1) Obj.XData(end)];
            yd=[Obj.YData(1) Obj.YData(end)];
            xstep=(xd(end)-xd(1))/size(C,2);
            ystep=(yd(end)-yd(1))/size(C,1);
            xd=xd+[-xstep/2 xstep/2];
            yd=yd+[-ystep/2 ystep/2];
            Obj.XData = linspace(xd(1),xd(2),size(C,2)+1);
            Obj.YData = linspace(yd(1),yd(2),size(C,1)+1)';
            Obj.ZData = zeros([size(C,1) size(C,2)]+1);
            Obj.EdgeColor = 'none';
            Obj.FaceColor = 'texturemap';
            Obj.FaceAlpha= 1;
            [F,V,Cface,Cedge,E,T]=surf2FV(Obj);
            Obj.Faces=F;
            Obj.E=E;
            Obj.Vertices=V;
            Obj.FaceVertexCData=Cface;
            Obj.EdgeVertexCData=Cedge;
            Obj.TextureVertices=T;
            data=addmesh(data,loc_scene,Obj);
        case 'light'
            data=addlight(data,loc_scene,Obj);
        case 'text'
            data=addtext(data,loc_scene,Obj);
        otherwise
            warning('figure2xhtml:figurobj',[type ' not supported']);
    end
end

function add_function_paths()
% add all needed function paths
try
    functionname='figure2xhtml.m';
    functiondir=which(functionname);
    functiondir=functiondir(1:end-length(functionname));
    addpath([functiondir '/functions'])
    addpath([functiondir '/functions_xml'])
catch me
    disp(me.message);
end

function hout=axis_handle(hin)
% Get axis handle of figure-handle or current handle
if(isempty(hin)),
    hout=gca;
else
    type=get(hin,'Type');
    if(~strcmpi(type,'axes'))
        Ch=get(hin,'Children');
        for i=1:length(Ch)
            type=get(Ch(i),'Type');
            if(strcmpi(type,'axes')), hout=Ch(i); break; end
        end
    else
        hout=hin;
    end
end


