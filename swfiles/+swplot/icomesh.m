function TR = icomesh(nSub)
% creates mesh by subdividing icosahedron faces
%
% TR = SWPLOT.ICOMESH(nSub)
%
% The output is a triangulated surface of the unit sphere, containing
% 20*4^nSub triangular faces. The output can be plotted using the trimesh()
% function.
%
% Input:
%
% nSub      Number of subdivisions. Default is 0 for icosahedron mesh
%           output.
%
% Output:
%
% TR        TriRep class triangulation object for plotting with trimesh().
%

% Using code from:
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: May.2012

if nargin == 0
    nSub = 0;
end

% get the vertex coordinates of the icosahedron
t = (1+sqrt(5))/2; % golden ratio
x = [0 1 t];
s = [1 1 1; 1 1 -1; 1 -1 -1; 1 -1 1];
x = repmat(x,[4 1]).*s;
x = [x;circshift(x,[0 -1]);circshift(x,[0 -2])];
xL2 = sqrt(sum(x.^2,2));
x   = bsxfun(@rdivide,x,xL2);

% use the result of convhulln(x) to speed up the calculation
cHullx = [1 5 9 6 4 9 1 4 10 4 1 9 8 2 5 8 1 10 1 8 5 12 6 9 5 12 9 2 12 5 12 ...
    3 6 3 12 2 3 7 6 4 7 10 7 4 6 8 11 2 11 3 2 11 8 10 7 11 10 11 7 3];
cHullx = reshape(cHullx,3,[])';
% triangulate the points
% only compatible with new Matlab versions
%TR = triangulation(fliplr(convhulln(x)),x);
TR = TriRep(fliplr(cHullx),x); %#ok<DTRIREP>

% get the data structure
%F = TR.ConnectivityList;
%X = TR.Points;
F = TR.Triangulation;
X = TR.X;

% spherical subdivision
for ii = 1:nSub
    % subdivide the mesh
    V1 = (X(F(:,1),:)+X(F(:,2),:))/2;
    V2 = (X(F(:,2),:)+X(F(:,3),:))/2;
    V3 = (X(F(:,3),:)+X(F(:,1),:))/2;
    V  = [V1;V2;V3];
    
    % remove repeating vertices
    % setOrder='stable' ensures that identical results (in terms of face
    %  connectivity) will be obtained for meshes with same topology
    [V,~,idx] = unique(V,'rows','stable');
    
    % assign indices to the new triangle vertices
    Nx = size(X,1); % # of vertices
    Nf = size(F,1); % # of faces
    
    V1 = Nx + idx(1:Nf);
    V2 = Nx + idx((Nf+1):2*Nf);
    V3 = Nx + idx((2*Nf+1):3*Nf);
    
    % define new faces
    T1 = [F(:,1) V1 V3];
    T2 = [F(:,2) V2 V1];
    T3 = [F(:,3) V3 V2];
    T4 = [V1     V2 V3];
    
    T1 = permute(T1,[3 1 2]);
    T2 = permute(T2,[3 1 2]);
    T3 = permute(T3,[3 1 2]);
    T4 = permute(T4,[3 1 2]);
    
    F = [T1;T2;T3;T4];
    F = reshape(F,[],3,1);
    
    % new mesh
    X  = cat(1,X,V);
    %TR = triangulation(F,X);
    %F  = TR.ConnectivityList;
    %X  = TR.Points;
    TR = TriRep(F,X); %#ok<DTRIREP>
    F  = TR.Triangulation;
    X  = TR.X;
    
end

if nSub > 0
    % reproject the points onto the surface of the unit sphere
    XL2 = sqrt(sum(X.^2,2));
    X   = bsxfun(@rdivide,X,XL2);
    %TR  = triangulation(F,X);
    TR  = TriRep(F,X);  %#ok<DTRIREP>
end

end