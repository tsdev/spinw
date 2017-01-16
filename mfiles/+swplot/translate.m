function translate(mode, hFigure)
% translate objects on swplot figure
%
% SWPLOT.TRANSLATE(mode, {hFigure})
%
% The function translates the objects of an swplot, where the coordinate
% system is defined by the plane of the figure with horizontal x-axis,
% vertical y-axis and out-of-plane z-axis.
%
% Input:
%
% mode      Either a vector with three elements determining the translation 
%           value in the figure plane coordinate system, or 'auto' that
%           centers the object on the figure. Default is 'auto'.
% hFigure   Handle of the swplot figure window, optional.
%
% See also SWPLOT.ZOOM.
%

if nargin == 0
    mode = 'auto';
end

if nargin < 2
    % find active figure
    hFigure = swplot.activefigure;
end

if isnumeric(mode)    
    T = swplot.transform(hFigure);
    %T(1:3,4) = mode(:);
    T(1:3,4) = T(1:3,1:3)\mode(:);
    swplot.transform(T,hFigure);
elseif strcmpi(mode,'auto')
    % find outer positions of objects
    obj = getappdata(hFigure,'objects');
    % position in lu
    pos = cat(3,obj(:).position);
    center = mean([min(pos(:,1,:),[],3) max(pos(:,1,:),[],3)],2);
    % translate in xyz units
    if swplot.ishg(hFigure)
        T = swplot.transform(hFigure);
    else
        T = eye(4);
    end
    T(1:3,4) = -swplot.base(hFigure)*center;
    swplot.transform(T,hFigure);

else
    error('translate:WrongInput','Wrong zoom mode!');
end

end