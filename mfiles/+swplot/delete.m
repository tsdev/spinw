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

% remove facepatch and edgepatch from handle list
fPatch = getappdata(hFigure,'facepatch');
ePatch = getappdata(hFigure,'edgepatch');

% all handles to delete
handle = [sObj(pIdx).handle];

% empty face patch if necessary
if ismember(fPatch,handle)
    fn = getappdata(fPatch,'facenumber');
    F  = get(fPatch,'Faces');
    V  = get(fPatch,'Vertices');
    C  = get(fPatch,'FaceVertexCData');
    % delete faces and vertices
    fToDel = ismember(fn,number);
    vToDel = F(fToDel,:);
    F(fToDel,:) = [];
    C(fToDel,:) = [];
    V(vToDel(:),:) = [];
    
    % renumber the face indices change old values to new values
    newF = changem(F,1:size(V,1),unique(F(:)'));
    set(fPatch,'Faces',newF,'Vertices',V,'FaceVertexCData',C);
    setappdata(fPatch,'facenumber',fn(~fToDel));
    % remove from handle list
    handle(handle==fPatch) = [];
end

% empty edge patch if necessary
if ismember(ePatch,handle)
    vn = getappdata(ePatch,'vertexnumber');
    F  = get(ePatch,'Faces');
    V  = get(ePatch,'Vertices');
    C  = get(ePatch,'FaceVertexCData');
    % delete faces and vertices
    vToDel = ismember(vn,number);
    fToDel = ismember(F(:,1),find(vToDel));
    F(fToDel,:) = [];
    C(vToDel,:) = [];
    V(vToDel,:) = [];
    % renumber the face indices change old values to new values
    newF = changem(F,1:size(V,1),unique(F(:)'));
    set(ePatch,'Faces',newF,'Vertices',V,'FaceVertexCData',C);
    setappdata(ePatch,'vertexnumber',vn(~vToDel));
    % remove from handle list
    handle(handle==ePatch) = [];
end

% delete objects
delete(handle);
setappdata(hFigure,'objects',sObj(~pIdx));

end