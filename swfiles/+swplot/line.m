function hPatch = line(varargin)
% creates 3D line patch
% 
% ### Syntax
% 
% `hLine = swplot.line(rStart, rEnd)`
%
% `hLine = swplot.line(r,[])`
% 
% `hLine = swplot.line(rStart, rEnd, lineStyle, lineWidth, multiPatch)`
%
% ### Description
%
% `hLine = swplot.line(rStart, rEnd)` creates disconnected line segments
% between multiple `rStart(:,i)` `rEnd(:,i)` pairs of 3D coordinates. The
% lines are shown as patch faces.
%  
% `hLine = swplot.line(r,[])` creates connected line segments  between
% the consicutive points `r(:,i)`.
%  
% `hPatch = swplot.line(handle, ...)` adds the generated patch object to a
% given axis if `handle` is an axis handle or adds the lines to an
% existing [matlab.patch] object, if the given `handle` points to a patch
% object.  
% 
% ### Input Arguments
% 
% `handle`
% : Handle of an axis or triangulated patch object. In case of patch
%   object, the constructed faces will be added to the existing object.
% 
% `rStart`
% : Coordinate(s) of the starting point, either a 3 element vector or
%   a matrix with dimensions of $[3\times n_{lineSegment}] to plot multiple line
%   segments.
% 
% `rEnd`
% : Coordinate(s) of the end point, either a 3 element vector or
%   a matrix with dimensions of $[3\times n_{lineSegment}]$ to plot multiple line
%   segments.
% 
% `r`
% : Matrix with dimensions of $[3\times n_{obj}\times n_{lineSegment}]$. The function
%   will plot $n_{obj}$ number of disconnected curves. The $i$th
%   curve will follow the `x=r(1,i,:)`, `y=r(2,i,:)`, `z=r(3,i,:)`
%   (parameteric) segmented curve.
% 
% `lineStyle`
% : Line style, default value is `'-'` for continuous line. Any other
%   Matlab line style string is accepted: `'--'`\|`':'`\|`'-.'`\|`'none'`.
% 
% `lineWidth`
% : Line width in pt, default value is 0.5.
% 
% `mPatch`
% : If `true`, a separate patch object will be created per line
%   segment. Default is `false`, a single patch object will store all
%   line segments.
% 
% ### See Also
% 
% [matlab.line]
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
    nArgExt = nargin-2;
    argExt  = {varargin{3:end}};
    
else
    hAxis   = gca;
    hPatch  = [];
    rStart  = varargin{1};
    
    nArgExt = nargin-1;
    argExt  = {varargin{2:end}};
    
end

if nArgExt > 0
    rEnd    = argExt{1};
else
    rEnd = [];
end

if nArgExt > 1
    lineStyle = argExt{2};
end
if nArgExt > 2
    lineWidth = argExt{3};
end
if nArgExt > 3
    mPatch = argExt{4};
end

if numel(rStart) == 3
    rStart = rStart(:);
    rEnd   = rEnd(:);
end
    
if ~isempty(rEnd)
    r = cat(3,rStart,rEnd);
else
    r = rStart;
end
    
if size(r,1)~=3
    error('line:WrongInput','To plot multiple line segments use matrices with dimensions of [3 nSegment]!');
end

% number of line segments
nObject = size(r,2);
% number of vertices per curve
nPoint = size(r,3);

% create the vertices
V = reshape(permute(r,[3 2 1]),[],3);
%V = reshape(permute(cat(3,rStart,rEnd),[1 3 2]),3,[])';

L = bsxfun(@plus,(1:(nPoint-1))',(0:(nObject-1))*nPoint);
F = [L(:) L(:)+1];

% black color
C = repmat([0 0 0],[size(V,1) 1]);
% transparency
A = ones(size(V,1),1);

if isnumeric(lineStyle) || numel(lineWidth)>1
    % multiple patch forced
    mPatch = true;
end

% all allowed linse style
lineStyle0 = {'-' '--' '-.' ':' 'none'};
% line style string in a cell
if ischar(lineStyle)
    lineStyle = {lineStyle};
else
    lineStyle = lineStyle0(lineStyle);
end

if isempty(hPatch)
    % create new patch
    if mPatch
        % prepare lineStyle
        if numel(lineStyle) == 1
            lineStyle = repmat(lineStyle,[1 nObject]);
        end
        if numel(lineWidth) == 1
            lineWidth = repmat(lineWidth,[1 nObject]);
        end

        % there is no zero linewidth
        lineWidth(lineWidth<=0) = 1e-3;
        
        hPatch = gobjects(1,nObject);
        for ii = 1:nObject
            hPatch(ii) = patch('Parent',hAxis,'Vertices',V((1:nPoint)+(ii-1)*nPoint,:),...
                'Faces',F((1:(nPoint-1)),:),'FaceLighting','none','AlphaDataMapping','none',...
                'EdgeColor',C(ii,:),'FaceColor','none','Tag','line',...
                'LineStyle',lineStyle{ii},'LineWidth',lineWidth(ii),...
                'EdgeAlpha',A(ii));
        end

    else
    hPatch = patch('Parent',hAxis,'Vertices',V,'Faces',F,'FaceLighting','none',...
        'EdgeColor','flat','FaceColor','none','Tag','line','AlphaDataMapping','none',...
        'FaceVertexCData',C,'LineStyle',lineStyle{1},'LineWidth',lineWidth,...
        'EdgeAlpha','flat','FaceVertexAlphaData',A);
    end
else
    % add to existing patch
    V0 = get(hPatch,'Vertices');
    F0 = get(hPatch,'Faces');
    C0 = get(hPatch,'FaceVertexCData');
    A0 = get(hPatch,'FaceVertexAlphaData');
    % number of existing faces
    nV0 = size(V0,1);
    set(hPatch,'Vertices',[V0;V],'Faces',[F0;F+nV0],'FaceVertexCData',[C0;C],...
        'FaceVertexAlphaData',[A0;A]);
end

if strcmp(get(hAxis,'Tag'),'swaxis') && ~mPatch
    % replicate the arrow handle to give the right number of added objects
    hPatch = repmat(hPatch,[1 nObject]);
end

end