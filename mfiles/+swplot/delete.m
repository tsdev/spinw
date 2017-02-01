function delete(varargin)
% deletes objectsand their data on swplot figure
%
% swplot.delete(hFigure,number)
%
% Deletes objects and their data that corresponds to the given unique
% numbers on hFigure swplot figure.
%
% swplot.delete(number)
%
% Deletes objects and their data that corresponds to the given unique
% numbers on the active swplot figure.
%
% If number equals to 0, all objects will be deleted from the figure.
%

if nargin == 1
    hFigure = swplot.activefigure;
    number  = varargin{1};
elseif nargin == 2
    hFigure = varargin{1};
    number  = varargin{2};
end

% get the objects
sObj = getappdata(hFigure,'objects');

if any(number==0)
    % select all objects to delete
    pIdx = 1:numel(sObj);
    number = 1:max([sObj(:).number]);
    % delete legend
    swplot.legend('off',hFigure)
    setappdata(hFigure,'legend',struct('handle',gobjects(0),'text',{''},...
        'type',[],'color',[],'name',{''}));
    
else
    % find objects with the given numbers
    pIdx = ismember([sObj(:).number],number);
end

% all handles to delete
handle = unique([sObj(pIdx).handle]);

% empty face patches if necessary
for ii = 1:numel(handle)
    hPatch = handle(ii);
    if isappdata(hPatch,'facenumber')
        fn = getappdata(hPatch,'facenumber');
        F  = get(hPatch,'Faces');
        V  = get(hPatch,'Vertices');
        C  = get(hPatch,'FaceVertexCData');
        A  = get(hPatch,'FaceVertexAlphaData');
        % delete faces and vertices
        fToDel = ismember(fn,number);
        
        if all(fToDel)
            delete(hPatch);
        elseif any(fToDel)
            vToDel = F(fToDel,:);
            F(fToDel,:) = [];
            C(fToDel,:) = [];
            A(fToDel,:) = [];
            V(vToDel(:),:) = [];
            % renumber the face indices change old values to new values
            newF = changem(F,1:size(V,1),unique(F(:)'));
            set(hPatch,'Faces',newF,'Vertices',V,'FaceVertexCData',C,'FaceVertexAlphaData',A);
            setappdata(hPatch,'facenumber',fn(~fToDel));
        end
    elseif isappdata(hPatch,'vertexnumber')
        vn = getappdata(hPatch,'vertexnumber');
        F  = get(hPatch,'Faces');
        V  = get(hPatch,'Vertices');
        C  = get(hPatch,'FaceVertexCData');
        A  = get(hPatch,'FaceVertexAlphaData');
        % delete faces and vertices
        vToDel = ismember(vn,number);
        if all(vToDel)
            delete(hPatch);
        elseif any(vToDel)
            fToDel = ismember(F(:,1),find(vToDel));
            F(fToDel,:) = [];
            C(vToDel,:) = [];
            A(vToDel,:) = [];
            V(vToDel,:) = [];
            % renumber the face indices change old values to new values
            newF = changem(F,1:size(V,1),unique(F(:)'));
            set(hPatch,'Faces',newF,'Vertices',V,'FaceVertexCData',C,'FaceVertexAlphaData',A);
            setappdata(hPatch,'vertexnumber',vn(~vToDel));
        end
    else
        delete(hPatch);
    end
end

setappdata(hFigure,'objects',sObj(~pIdx));

end