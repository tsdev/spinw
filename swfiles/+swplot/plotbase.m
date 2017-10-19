function varargout = plotbase(varargin)
% plots basis vectors
% 
% ### Syntax
% 
% `swplot.plotbase(Name,Value)`
% 
% `hFigure = swplot.plotbase(Name,Value)`
%
% ### Description
% 
% `swplot.plotbase(Name,Value)` plots the three basis vectors that define
% the coordinate system of the plot, either the $abc$ lattice vectors,
% $xyz$ Descartes coodinate system or the $hkl$ reciprocal lattice vectors.
% 
% ### Name-Value Pair Arguments
% 
% `'mode'`
% : String that determines the type of basis vectors to plot. Possible
%   values are:
%   * `abc`     plots the lattice vectors (default),
%   * `hkl`     plots the reciprocal lattice vectors,
%   * `xyz`     plots the $xyz$ Descartes coordinate system.
% 
% `'length'`
% : Determines the length of the 3 basis vectors. If 0, the
%   length won't be rescaled. If non-zero, the `length` parameter
%   determines the length of the plotted vectors in \\ang. Default value is
%   2 \\ang.
% 
% `'label'`
% : Logical variable, plots the vector labels if `true`. Default value is 
%   `true`.
% 
% `'color'`
% : Color of the arrows, either a cell of three color name strings or a
%   matrix with dimensions of ${3\times 3]$ where each column defines the
%   RGB values of a color. Default value is `{'red' 'green' 'blue'}`.
% 
% `'R'`
% : Radius of the arrow body, default value is 0.06 \\ang.
% 
% `'alpha'`
% : Head angle of the arrow in degree units, default value is 30\\deg.
% 
% `'lHead'`
% : Length of the arrow head, default value is 0.5 \\ang.
% 
% `'d'`
% : Distance from origin in $xyz$ units, default value is `[1 1 1]`.
% 
% `'dtext'` : Distance of the label from the arrow in xyz units, default
%   value is 0.5 \\ang.
% 
% #### General paraters
%
% These parameters have the same effect on any of the `swplot.plot...`
% functions.
%
% `'obj'`
% : [spinw] object.
% 
% `'figure'`
% : Handle of the swplot figure. Default value is the active figure handle.
% 
% `'tooltip'`
% : If `true`, the tooltips will be shown when clicking on the plot
%   objects. Default value is `true`.
% 
% `'shift'`
% : Column vector with 3 elements, all object positions will be
%   shifted by the given value in \\ang units. Default value is
%   `[0;0;0]`.
% 
% `'replace'`
% : If `true` the plot will replace the previous plot of the same type.
%   Default value is `true`.
% 
% `'translate'`
% : If `true`, the plot will be centered, independent of the range. Default
%   value is `false`.
% 
% `'zoom'`
% : If `true`, the swplot figure will be zoomed to make the plot objects
%   cover the full figure. Default is `true`.
% 
% `'copy'`
% : If `true`, a clone of the [spinw] object will be saved in the
%   swplot figure data which can be retwrived using
%   `swplot.getdata('obj')`. If `false`, the handle of the original [spinw]
%   object is saved which is linked to the input `obj` and so it changes
%   when `obj` is changed. Default value is `false`.
%
% `nPatch`
% : Number of vertices on any patch object that is not the icosahedron,
%   default value is stored in `swpref.getpref('npatch')`.
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