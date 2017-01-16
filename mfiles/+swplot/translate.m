function translate(mode, hFigure)
% translate objects on swplot figure
%
% SWPLOT.TRANSLATE(mode, {hFigure})
%
% Input:
%
% mode      Either a vector with three elements determining the translation 
%           value in xyz units, or 'auto' that centers the object on the
%           figure. Default is 'auto'.
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
    T(1:3,4) = mode(:);
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