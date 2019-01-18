function hPatch = polyhedron(varargin)
% creates convex polyhedra or polygon from vertex list
% 
% ### Syntax
% 
% `hPatch = swplot.polyhedron(vertices)`
% 
% `hPatch = swplot.polyhedron(handle, ...)`
%
% ### Description
% 
% `hPatch = swplot.polyhedron(vertices)` creates convex polyhedra or
% polygon from a given vertex list (unordered list of 3D coordinates). To
% draw the polyhedron the convex hull of the given point cloud is
% calculated using [matlab.convhulln]. It is automatically detected if the
% given vertex points lie on a plane in which case the convex polygon is
% drawn.
%  
% `hPatch = swplot.polyhedron(handle, ...)` adds the generated patch object
% to a given axis if `handle` is an axis handle or adds the polyhedron to
% an existing [matlab.patch] object, if the given `handle` points to a
% patch object.
% 
% ### Input Arguments
% 
% `vertices`
% : Matrix with dimensions of $[3\times n_{obj}\times n_{point}]$, where
%   $n_{obj}$ is the number of polyhedra to draw, $n_{point}$ is the number
%   of vertices per polyhedron.
% 
% ### See Also
% 
% [matlab.convhulln]
%

if nargin == 0
    swhelp swplot.polyhedron
    return
end

if numel(varargin{1}) == 1
    % first input figure/patch handle
    if strcmp(get(varargin{1},'Type'),'axes')
        hAxis  = varargin{1};
        hPatch = [];
    else
        hAxis  = gca;
        hPatch = varargin{1};
    end
    
    vertices  = varargin{2};
    
else
    hAxis  = gca;
    hPatch = [];
    vertices  = varargin{1};
end

vertices = permute(vertices,[1 3 2]);

if size(vertices,1) ~=3 || size(vertices,2)<3
    error('polyhedron:WrongInput','The vertex matrix dimension has to be [3 nObject nPoint>=3]!');
end

nObject = size(vertices,3);
nPoint  = size(vertices,2);

F = cell(1,nObject);

for ii = 1:nObject
    V = vertices(:,:,ii);
    switch rank(bsxfun(@minus,V,V(:,1)))
        case 3
            % calculate the faces of the 3D polyhedra
            F{ii}  = convhulln(V');
            
        case 2
            % flat polygons
            % normal to the plane of the polygon
            base = orth(bsxfun(@minus,V,V(:,1)));
            v1 = base(:,1);
            v2 = base(:,2);
            v3 = cross(v1,v2);
            
            T   = [v1 v2 v3];
            V2D = T'*V;
            F2D = convhulln(V2D(1:2,:)');
            F{ii}   = [F2D(1)*ones(size(F2D,1)-2,1) F2D(2:(end-1),1) F2D(3:(end),1)];
        case 1
            error('polyhedron:WrongInput','The given vertices are along a line!');
        case 0
            error('polyhedron:WrongInput','The given vertices are identical!')
    end
end

% check if a single patch object can be used
nF = cellfun(@(C)numel(C),F);

if any(nF-nF(1))
    mPatch = true;
else
    mPatch = false;
end

if mPatch
    % plot each polyhedron in separate patch
    hPatch = gobjects(1,nObject);
    for ii = 1:nObject
        % create the patch
        F0 = F{ii};
        nF = size(F0,1);
        hPatch(ii) = patch('Parent',hAxis,'Vertices',vertices(:,:,ii)','Faces',F0,...
            'FaceLighting','flat','EdgeColor','none','FaceColor','flat',...
            'Tag','polyhedron','AlphaDataMapping','none','FaceAlpha','flat',...
            'FaceVertexAlphaData',ones(nF,1),'FaceVertexCData',repmat([1 0 0],[nF 1]));
        
    end
else
    % create the face and vertex list for a single patch
    V = reshape(vertices,3,[])';
    F = cell2mat(permute(F,[1 3 2]));
    F = bsxfun(@plus,F,permute((0:(nObject-1))*nPoint,[1 3 2]));
    F = reshape(permute(F,[2 1 3]),3,[])';
    nF = size(F,1);
    A = ones(nF,1);
    C = repmat([1 0 0],[nF 1]);
    
    if isempty(hPatch)
        hPatch = patch('Parent',hAxis,'Vertices',V,'Faces',F,'FaceLighting','flat',...
            'EdgeColor','none','FaceColor','flat','Tag','polyhedron','AlphaDataMapping','none',...
            'FaceAlpha','flat','FaceVertexAlphaData',A,'FaceVertexCData',C);
    else
        V0 = get(hPatch,'Vertices');
        F0 = get(hPatch,'Faces');
        C0 = get(hPatch,'FaceVertexCData');
        A0 = get(hPatch,'FaceVertexAlphaData');
        % number of existing faces
        nV0 = size(V0,1);
        set(hPatch,'Vertices',[V0;V],'Faces',[F0;F+nV0],'FaceVertexCData',[C0;C],...
            'FaceVertexAlphaData',[A0;A]);
    end
    
    if strcmp(get(hAxis,'Tag'),'swaxis')
        % replicate the arrow handle to give the right number of added objects
        % for swplot figure
        hPatch = repmat(hPatch,[1 nObject]);
    end
end

end