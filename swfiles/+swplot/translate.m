function translate(mode, hFigure)
% translates objects on swplot figure
% 
% ### Syntax
% 
% `swplot.translate(mode)`
% 
% `swplot.translate(mode, hFigure)`
%
% ### Description
% 
% `swplot.translate(mode)` translates the objects of an active swplot
% figure, where the coordinate system is defined by the plane of the figure
% with horizontal $x$-axis, vertical $y$-axis and out-of-plane $z$-axis.
% 
% `swplot.translate(mode, hFigure)` acts on the figure referenced by the
% `hFigure` handle.
%
% ### Input Arguments
% 
% `mode`
% : Either a vector with three numbers that determine the translation 
%   vector in the figure plane coordinate system, or `'auto'` that
%   centers figure to the middle of the objects. Default value is `'auto'`.
% 
% `hFigure`
% : Handle of the swplot figure, default value is the active figure.
% 
% ### See Also
% 
% [swplot.zoom]
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
    if any(isnan(center))
        return
    end
    
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