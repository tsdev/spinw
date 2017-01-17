function hPatch = line(varargin)
% draws a 3D line using patch
%
% hLine = SWPLOT.LINE(rStart, rEnd, {lineStyle}, {lineWidth},{multiPatch})
%
% hLine = SWPLOT.LINE(handle,...)
%
% Handle can be the handle of an axes object or a line object. It either
% selects an axis to plot or a patch object (triangulated) to add vertices
% and faces.
%
% Input:
%
% handle    Handle of an axis or patch object. In case of patch object, the
%           constructed faces will be added to the existing object instead
%           of creating a new one.
% rStart    Coordinate(s) of the starting point, either a 3 element vector or
%           a matrix with dimensions [3 nLineSegment] to plot multiple line
%           segments.
% rEnd      Coordinate(s) of the end point, either a 3 element vector or
%           a matrix with dimensions [3 nLineSegment] to plot multiple line
%           segments.
% lineStyle Line style, default is '-' for continuous line.
% lineWidth Line with in pt, default is 0.5.
% mPatch    If true, a separate patch object will be created per line
%           segment. Default is false, a single patch object will store all
%           line segments.
%
% See also LINE.
%

if nargin == 0
    help swplot.line
    return
end

lineStyle = '-';
lineWidth = 0.5;
mPatch    = false;

if numel(varargin{1}) == 1
    % first input figure/patch handle
    if strcmp(get(varargin{1},'Type'),'axes')
        hAxis  = varargin{1};
        hPatch = [];
    else
        hAxis  = gca;
        hPatch = varargin{1};
    end
    rStart  = varargin{2};
    rEnd    = varargin{3};
    nArgExt = nargin-3;
    argExt  = {varargin{4:end}};
    
else
    hAxis   = gca;
    hPatch  = [];
    rStart  = varargin{1};
    rEnd    = varargin{2};
    nArgExt = nargin-2;
    argExt  = {varargin{3:end}};
    
end

if nArgExt > 0
    lineStyle = argExt{1};
end
if nArgExt > 1
    lineWidth = argExt{2};
end
if nArgExt > 2
    mPatch = argExt{3};
end

if numel(rStart) == 3
    rStart = rStart(:);
    rEnd   = rEnd(:);
else
    if size(rStart,1)~=3 || size(rEnd,1)~=3
        error('line:WrongInput','To plot multiple line segment use [3 nSegment] matrices!');
    end
end

% number of line segments
nObject = size(rStart,2);

% create the vertices
V = reshape(permute(cat(3,rStart,rEnd),[1 3 2]),3,[])';

L = (1:nObject)';
F = [2*L-1 2*L];

% black color
C = repmat([0 0 0],[size(V,1) 1]);

if isnumeric(lineStyle) || numel(lineWidth)>1
    % multiple patch forced
    mPatch = true;
end

if strcmp(get(hAxis,'Tag'),'swaxis') && ~mPatch && strcmp(lineStyle,'-') && lineWidth == 0.5
    % make sure we are on the plot axis of an swobject
    % add object to the existing edge patch, but only for continuous lines
    hFigure = get(hAxis,'Parent');
    hPatch = getappdata(hFigure,'edgepatch');
end

if isempty(hPatch)
    % create new patch
    if mPatch
        % prepare lineStyle
        if ischar(lineStyle)
            lineStyle = numel(lineStyle);
        end
        if numel(lineStyle) == 1
            lineStyle = repmat(lineStyle,[1 nObject]);
        end
        if numel(lineWidth) == 1
            lineWidth = repmat(lineWidth,[1 nObject]);
        end

        hPatch = gobjects(1,nObject);
        for ii = 1:nObject
            hPatch(ii) = patch('Parent',hAxis,'Vertices',V(2*ii+[-1 0],:),...
                'Faces',[1 2],'FaceLighting','none',...
                'EdgeColor',C(ii,:),'FaceColor','none','Tag','line',...
                'LineStyle',repmat('-',[1 lineStyle(ii)]),'LineWidth',lineWidth(ii));
        end

    else
    hPatch = patch('Parent',hAxis,'Vertices',V,'Faces',F,'FaceLighting','none',...
        'EdgeColor','flat','FaceColor','none','Tag','line',...
        'FaceVertexCData',C,'LineStyle',lineStyle,'LineWidth',lineWidth);
    end
else
    % add to existing patch
    V0 = get(hPatch,'Vertices');
    F0 = get(hPatch,'Faces');
    C0 = get(hPatch,'FaceVertexCData');
    % number of existing faces
    nV0 = size(V0,1);
    set(hPatch,'Vertices',[V0;V],'Faces',[F0;F+nV0],'FaceVertexCData',[C0;C]);
end

if strcmp(get(hAxis,'Tag'),'swaxis') && strcmp(lineStyle,'-') && lineWidth == 0.5 && ~mPatch
    % replicate the arrow handle to give the right number of added objects
    hPatch = repmat(hPatch,[1 nObject]);
end

end