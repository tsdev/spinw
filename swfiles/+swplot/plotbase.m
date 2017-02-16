function varargout = plotbase(varargin)
% plots the edges of unit cells on swplot figure
%
% SWPLOT.PLOTBASE('Option1',Value1,...)
%
% hFigure = SWPLOT.PLOTBASE('Option1',Value1,...)
%
% Options:
%
% mode      String that determines the type base to plot. Possible values
%           are:
%               abc     Plots the lattice vectors, default.
%               xyz     Plots the lattice Descartes coordinate system.
%               hkl     Plots the reciprocal lattice vectors.
% length    Determines the length of the a, b and c arrows. If 0, the
%           length will be equal to the corresponding lattice parameters,
%           while if non-zero, the number determines the length in
%           Angstrom. Default is 2 Angstrom.
% label     Logical variable, plots abc labels if true. Default is true.
% figure    Handle of the swplot figure. Default is the selected figure.
% color     Color of the arrows, default is red-green-blue for abc, stored
%           in the columns of a 3x3 matrix.
% R         Radius value of arrow body, default is 0.06.
% alpha     Head angle of the arrow in degree units, default is 30 degree.
% lHead     Length of the arrow head, default value is 0.5.
% d         Distance from origin in xyz units.
% dtext     Distance from arrow and in xyz units.
% shift     Column vector with 3 elements, the basis vectors will be
%           shifted by the given values in Angstrom unit. Default value is
%           [0;0;0].
% translate If true, all plot objects will be translated to the figure
%           center. Default is false.
% zoom      If true, figure will be automatically zoomed to the ideal size.
%           Default is false.
%

% default values
col0    = swplot.color({'red' 'green' 'blue'});
d0      = ones(1,3);
nMesh0  = swpref.getpref('nmesh',[]);
nPatch0 = swpref.getpref('npatch',[]);

inpForm.fname  = {'range' 'mode'    'figure' 'color' 'R'   'alpha' 'lhead' 'shift'};
inpForm.defval = {[]      'abc'     []       col0    0.06  30      0.5     [0;0;0]};
inpForm.size   = {[-1 -2] [1 -3]    [1 1]    [-4 -5] [1 1] [1 1]   [1 1]   [3 1]  };
inpForm.soft   = {true    false     true     false   false false   false   false  };

inpForm.fname  = [inpForm.fname  {'d'   'dtext' 'label' 'tooltip' 'translate' 'copy'}];
inpForm.defval = [inpForm.defval {d0    0.5     true    true      false       false }];
inpForm.size   = [inpForm.size   {[1 3] [1 1]   [1 1]   [1 1]     [1 1]       [1 1] }];
inpForm.soft   = [inpForm.soft   {false false   false   false     false       false }];

inpForm.fname  = [inpForm.fname  {'zoom' 'replace' 'nmesh' 'npatch' 'unit' 'length' 'legend'}];
inpForm.defval = [inpForm.defval {false  true      nMesh0  nPatch0  'lu'   2        true    }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]     [1 1]   [1 1]    [1 -6] [1 1]    [1 1]   }];
inpForm.soft   = [inpForm.soft   {false  false     false   false    false  false    false   }];

param = sw_readparam(inpForm, varargin{:});

if isempty(param.figure)
    hFigure  = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

% select range
if numel(param.range) == 6
    range = param.range;
elseif isempty(param.range)
    % get range from figure
    fRange = getappdata(hFigure,'range');
    if isempty(fRange)
        % fallback to default range
        range = [0 1;0 1;0 1];
    else
        % get plotting range and unit
        range       = fRange.range;
        param.unit  = fRange.unit;
    end
elseif numel(param.range) == 3
    % change range, if the number of unit cells are given
    range = [zeros(3,1) param.range(:)];
else
    error('plotbase:WrongInput','The given plotting range is invalid!');
end

% lower range limit, where basis vectors come
range0 = range(:,1);

% basis vectors
BV = swplot.base(hFigure);

switch param.mode
    case 'abc'
        pBase = eye(3);
        axText = {'a' 'b' 'c'};
    case 'xyz'
        pBase = inv(BV);
        axText = {'x' 'y' 'z'};
    case 'hkl'
        pBase = 2*pi*inv(BV)*inv(BV)';
        axText = {'h' 'k' 'l'};
end

switch param.unit
    case 'lu'
        % do nothing
    case 'xyz'
        % convert to lu
        range0 = BV\range0;
    otherwise
        error('plotbase:WrongInput','Option unit has invalid value!');
end

% convert d from xyz to base
d = (param.d/(BV'))'-range0;

% vector length
length0 = sqrt(sum((BV*pBase).^2,1));

if param.length == 0
    % do nothing
else
    % normalize the length of the vectors
    pBase  = bsxfun(@times,pBase,param.length./length0);
    length0 = param.length;
end

pos = bsxfun(@minus,cat(3,zeros(3),pBase),d);

% shift
pos = bsxfun(@plus,pos,BV\param.shift);

% generate data
data = repmat({BV},[1 3]);

if param.label
    % convert d from xyz to base
    pos_text = bsxfun(@minus,bsxfun(@times,pBase,1+param.dtext./length0),d);
    pos_text = bsxfun(@plus,pos_text,BV\param.shift);
    
    swplot.plot('type','text','position',pos_text,'text',axText,'data',data,...
        'figure',hFigure,'legend',false,'tooltip',false,'translate',false,...
        'zoom',false,'name','base_label','replace',param.replace);
end

% plot the arrows
swplot.plot('type','arrow','position',pos,'figure',hFigure,'color',param.color,...
    'name','base','legend',false,'tooltip',false,'translate',param.translate,...
    'zoom',param.zoom,'data',data,'replace',param.replace,'npatch',param.npatch);

% save range
setappdata(hFigure,'range',struct('range',range,'unit',param.unit));

if param.tooltip
    swplot.tooltip('on',hFigure);
end

if nargout>0
    varargout{1} = hFigure;
end

end