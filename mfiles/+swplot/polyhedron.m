function hPatch = polyhedron(varargin)
% draw convex polyhedra or polynom from vertex list
%
% hPatch = SWPLOT.POLYHEDRON(vertices)
%
% hPatch = SWPLOT.POLYHEDRON(handle,...
%
% Handle can be the handle of an axes object or a patch object. It either
% selects an axis to plot or a patch object (triangulated) to add vertices
% and faces.
%
% Input:
%
% vertices      Matrix with dimensions [3 nObject nPoint], where nObject is
%               the number of polyhedra to draw, nPoint is the number of
%               vertices per polyhedron.
%
% See also CONVHULLN.
%

if nargin == 0
    help swplot.polyhedron
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
    switch rank(V)
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
            'Tag','ellipsoid','AlphaDataMapping','none','FaceAlpha','flat',...
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
            'EdgeColor','none','FaceColor','flat','Tag','ellipsoid','AlphaDataMapping','none',...
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