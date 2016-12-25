function add(hAdd, hFigure)
% adds a graphical object to a figure
%
% SWPLOT.ADD(hAdd, {hFigure})
%
% It adds a graphical object to the hgtransform object of the figure for
% continuous rotation with the mouse.
%
% Input:
%
% hAdd      Struct, that contains the handles of the graphical objects,
%           e.g. hAdd.objtype1 = [handle1 handle2 ...].
% hFigure   The handle of the figure (number in the figure title bar). The
%           default is the active swplot figure.
%
% See also SWPLOT, SWPLOT.FIGURE, HGTRANSFORM.
%

if nargin == 0
    help swplot.add
    return
end

if nargin < 2
    % find active figure
    hFigure = swplot.activefigure;
end

cva = get(gca,'CameraViewAngle');
hTransform = getappdata(hFigure,'h');

if isappdata(hFigure,'objects')
    sObject = getappdata(hFigure,'objects');
else
    sObject = struct;
end

% convert a simple list of handles to the required structure to store in
% swplot
if ~isstruct(hAdd)
    nObjAdd = numel(hAdd);
    sTemp1 = struct('handle',cell(1,nObjAdd),'text',cell(1,nObjAdd),'position',cell(1,nObjAdd));
    sTemp2.general = sTemp1;
    % convert vector of graphical objects to cell
    hAdd = num2cell(hAdd(:));
    % save handles
    [sTemp2.general(:).handle] = hAdd{:};
    % empty placeholders
    text1     = repmat({''},nObjAdd,1);
    [sTemp2.general(:).text] = text1{:};
    position1 = repmat({nan(3,2)},nObjAdd,1);
    [sTemp2.general(:).position] = position1{:};
    hAdd = sTemp2;
end

% add all the objects to the hgtransform object.
names = fieldnames(hAdd);
for ii = 1:numel(names)
    hSelect = [hAdd.(names{ii})(:).handle];
    set(hSelect,'Parent',hTransform);
    set(hSelect,'Clipping','Off');
end

% comb together the handles of the old and new graphical objects.
namesAdd = fieldnames(hAdd);
for ii = 1:length(namesAdd)
    if isfield(sObject,namesAdd{ii})
        sObject.(namesAdd{ii}) = [sObject.(namesAdd{ii}) hAdd.(namesAdd{ii})];
    else
        sObject.(namesAdd{ii}) = hAdd.(namesAdd{ii});
    end
end

% Shift the origin to center the plot.
if isappdata(hFigure,'param')
    param = getappdata(hFigure,'param');
else
    param = struct;
end

% center object if it is crystal
% TODO change to BV matrix
if isfield(param,'range') && isappdata(hFigure,'obj')
    range       = param.range;
    basisVector = getappdata(hFigure,'obj');
    basisVector = basisVector.basisvector;
    T           = makehgtform('translate',-sum(basisVector * sum(range,2)/2,2)');
    set(hTransform,'Matrix',get(hTransform,'Matrix')*T);
end

% Saves the object handles into the figure UserData property.
setappdata(hFigure,'objects',sObject);
setappdata(hFigure,'h',hTransform);
set(gca,'CameraViewAngle',cva);

end