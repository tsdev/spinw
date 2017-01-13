function varargout = transform(varargin)
% transform objects on swplot figure
%
% SWPLOT.TRANSFORM(M, {hFigure})
%
% Transforms the objects on the active swplot figure using the
% transformation matrix M.
%
% Input:
%
% M         Transformation matrix, either has dimensions of 4x4 following
%           the Matlab standard for hgtransform, or 3x4 following the SpinW
%           format for space group transformations. Setting it to 0,
%           returns to the original orientation of the objects (equivalent
%           to eye(4)).
% hFigure   Handle of the swplot figure window, optional.
%
%
% M = SWPLOT.TRANSFORM({hFigure})
%
% Returns the transformation matrix (4x4) of the active swplot figure.
%
% If the figure is created without the hgtransform object, the
% transformation matrix moves the camera.
%
% See also SWPLOT.FIGURE, HGTRANSFORM.
%

M = [];

if nargin == 0
    % find active figure
    hFigure = swplot.activefigure;
elseif nargin == 1 && numel(varargin{1})==1 && varargin{1}~=0
    hFigure = varargin{1};
elseif nargin == 1 && (numel(varargin{1})>1 || varargin{1}==0)
    M       = varargin{1};
    hFigure = swplot.activefigure;
elseif nargin == 2
    M       = varargin{1};
    hFigure = varargin{2};
else
    error('transform:WrongInput','Wrong input, check help swplot.transform!');
end

h = getappdata(hFigure,'h');

if numel(M)==1 && M==0
    % generate unit matrix for special input 0
    M = eye(4);
end

if ~isempty(M) && all(size(M)==[3 4])
    M = [M;[0 0 0 1]];
end

if ~isempty(M) && ~all(size(M)==[4 4])
    error('transform:WrongInput','M matrix has wrong dimensions!');
end

if isempty(h)
    % nohg
    hAxis = getappdata(hFigure,'axis');
    if ~isempty(M)
        % rotate camera
        R = M(1:3,1:3)';
        set(hAxis,'CameraPosition', R*[0 0 100]',...
            'CameraTarget',         [0 0 0],...
            'CameraUpVector',       R*[0 1 0]');
    else
        %v1 = get(hAxis,'CameraPosition');
        %v2 = get(hAxis,'CameraUpVector');
        %
        %v = [v1;v2];
        %v10 = [0 0 norm(v1)];
        %v20 = [0 norm(v2) 0];
        %v0 = [v10;v20];
        warning('transform:MissingFeature','This option is not implemented yet!')
    end
else
    % hg
    if isempty(M)
        MT = get(h,'Matrix');
        MR = get(get(h,'Parent'),'Matrix');
        M = eye(4);
        M(1:3,1:3) = MR(1:3,1:3);
        M(1:3,4)   = MT(1:3,4);
        varargout{1} = M;
    else
        MT = eye(4);
        MT(1:3,4) = M(1:3,4);
        MR = eye(4);
        MR(1:3,1:3) = M(1:3,1:3);
        set(h,'Matrix',MT);
        set(get(h,'Parent'),'Matrix',MR);
    end
end

end