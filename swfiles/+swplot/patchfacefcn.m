function varargout = patchfacefcn(obj,hit,fun,selection,dLim,hTransform)
% callback function for patch face selection
% 
% ### Syntax
% 
% `patchfacefcn(hPatch, hit, callbackFun, selection)`
% 
% `patchfacefcn(hPatch, hit, callbackFun, selection, dLim)`
%
% ### Description
% 
% `patchfacefcn(hPatch, hit, callbackFun, selection)` finds the index of the
% face in a patch object which was clicked on by the mouse. The function
% should be used as a callback function for the `ButtonDownFcn` event of
% the patch object and it will call a user defined function with the same
% arguments as the `ButtonDownFcn` call, plus adding an extra argument, the
% face index. Thus the user defined callback function needs to have the
% following header:
% ``` 
% callbackFun(hPatch,hit,faceIndex)
% ```
%  
% The function can detect if the mouse click was on a face or on an edge of
% the patch object.
%  
% `patchfacefcn(hPatch, hit, callbackFun, selection, dLim)` adds optional
% control on the selectivity whether the click was on a face of a patch
% object.
%
% ### Examples
% 
% The color of any face of the red triangulated icosahedron will be
% changed from red to green if clicked on with the mouse.
%
% ```
% >>mesh = swplot.icomesh(1)
% >>V = mesh.X
% >>F = mesh.Triangulation
% >>hPatch = patch('Faces',F,'Vertices',V,'FaceColor','r','EdgeColor','none')
% >>axis equal
% >>axis off
% >>box on
% >>view(3)
% >>camlight('right')
% >>hPatch.FaceColor = 'flat'
% >>hPatch.FaceVertexCData = repmat([1 0 0],[size(F,1) 1])
% >>fun = @(hPatch,hif,face)set(hPatch,'FaceVertexCData',[hPatch.FaceVertexCData(1:(face-1),:); [0 1 0]; hPatch.FaceVertexCData((face+1):end,:)])
% >>hPatch.ButtonDownFcn = @(hPatch,hit)swplot.patchfacefcn(hPatch,hit,fun,'face')
% ```
% 
% ### Input Arguments
% 
% `hPatch`
% : Handle of the patch object.
% 
% `hit`
% : Hit object that contains the point where the object was hit.
% 
% `callbackFun`
% : User function that will be called in case of a click event
%     on `hPatch` object. It should have the following header:
%         `callbackFun(hPatch,hit,faceIndex)`
%     where `faceIndex` contains the index of the face that was
%     clicked on, it can contain a single index or more depending
%     on the selection type.
% 
% `selection`
% : String, that defines three diferent selection criteria when the
%   `callbackfun()` function will be called:
%   * `'face'`  The `callbackFun()` will be triggered if a face
%               was clicked (`numel(faceIndex)==1`).
%   * `'edge'`  The `callbackfun()` will be triggered if an edge
%               is clicked (`numel(faceIndex)==2`).
%   * `'all'`   The `callbackfun()` will be triggered for both
%               faces and edges.
% 
% `dLim`
% : Upper limit of the absolute value of the determinant that
%   determines whether a point is on a plane spanned by the two
%   edges of a triangle. Default value is $10^{-7}$ (tested).
% 
% ### See Also
% 
% [matlab.patch]
%

if isempty(hit)
    % older versions of Matlab hit is empty
    if nargout > 0
        varargout{1} = [];
    end
    return
end

% point where we hit the surface
P = hit.IntersectionPoint;

nF = size(obj.Faces,1);

% vertices of patch
V = obj.Vertices;

if size(V,2) ==2
    flat = true;
    % for 2D patch
    P = P(1:2);
else
    flat = false;
end

if size(obj.Faces,2) ~=3
    error('patchfacefcn:WrongObject','The pathcfacefcn() callback function works only for triangulated surfaces!');
end

% precision for finding planes of faces
if nargin < 5 || isempty(dLim)
    dLim = 1e-7;
    tLim = 0;
elseif numel(dLim==2)
    tLim = dLim(2);
    dLim = dLim(1);
else
    tLim = 0;
end

% limit on the edge of the triangle

if nargin<6 || isempty(hTransform)
    % no hgtransform parent, use unit transformation matrix
    M = eye(4);
else
    % get the transformation matrix
    M1 = get(hTransform,'Matrix');
    hTransform2 = get(hTransform,'Parent');
    M2 = get(hTransform2,'Matrix');
    %M = M1*M2;
    M = M2*M1;
end

% transform the point to the coordinate system of the object
T = M(1:3,4)';
R = M(1:3,1:3);
P = (P-T)*R;
if ~verLessThan('matlab','9.4')
    P = (P-T)*R;
end

% shift the origin to the first vertex of every triangle
E1 = V(obj.Faces(:,2),:)-V(obj.Faces(:,1),:);
E2 = V(obj.Faces(:,3),:)-V(obj.Faces(:,1),:);
D  = bsxfun(@minus,P,V(obj.Faces(:,1),:));

if ~flat
    % check whether the plane of any face contains the point
    det = sum(cross(D,E1,2).*E2,2);
    
    % find points within the plane
    pIdx = find(abs(det)<dLim);
    
    % determine barycentric coordinates
    bCoord = zeros(numel(pIdx),2);
    for ii = 1:numel(pIdx)
        bCoord(ii,:) = ([E1(pIdx(ii),:)' E2(pIdx(ii),:)']\D(pIdx(ii),:)')';
    end
    
    % find the right face(s)
    fIdx = pIdx(all(bCoord>=-tLim & bCoord<=1+tLim & repmat(sum(bCoord,2)<=1+tLim,1,2),2));
else
    % flat patch, all triangles are in the plane
    % determine barycentric coordinates
    bCoord = zeros(nF,2);
    for ii = 1:nF
        bCoord(ii,:) = ([E1(ii,:)' E2(ii,:)']\D(ii,:)')';
    end
    
    % find the right face(s)
    fIdx = find(all(bCoord>=0 & bCoord<=1 & repmat(sum(bCoord,2)<=1,1,2),2));
    
end

% % for debugging
% if isempty(fIdx)
%     error
% end

if isempty(fun)
    % don't call a function, just return the index of the patch face
    varargout{1} = fIdx;
else
    % call user defined callback function
    switch selection
        case 'all'
            % include clicks on edges (numel(fIdx)>1)
            if numel(fIdx) > 0
                fun(obj,hit,fIdx);
            end
        case 'face'
            % only trigger for faces
            if numel(fIdx) == 1
                fun(obj,hit,fIdx);
            end
        case 'edge'
            % only trigger for faces
            if numel(fIdx) == 2
                fun(obj,hit,fIdx);
            end
        otherwise
            error('pathcfacefcn:WrongInput','The pathcfacefcn() callback has only two modes: ''body'' and ''face''!');
    end
end

end