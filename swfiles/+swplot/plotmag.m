function varargout = plotmag(varargin)
% plots magnetic structure
% 
% ### Syntax
% 
% `swplot.plotmag(Name,Value)`
% 
% `hFigure = swplot.plotmag(Name,Value)`
% ### Description
% 
% `swplot.plotmag(Name,Value)` plots the magnetic structure stored in a
% [spinw] object onto an swplot figure. The magnetic structure is
% represented by arrows on the magnetic atoms.
% 
% ### Name-Value Pair Arguments
% 
% `'mode'`
% : String that defines the way the magnetic moments are plotted:
%   * `'all'`       Plot both the rotation plane of single-k incommensurate
%                   magnetic structures and the moment directions.
%   * `'circle'`    Plot only the rotation plane of incommensurate
%                   magnetic structures.
%   * `'arrow'`     Plots only the moment directions.
% 
% `'label'`
% : Whether to plot labels for magnetic atoms, default value is `true`.
% 
% `'dText'`
% : Distance between atom and its text label, default value is 0.1
%   \\ang.
% 
% `'fontSize'`
% : Font size of the text labels in pt, default value is stored in
%   `swpref.getpref('fontsize')`.
% 
% `'color'`
% : Color of the magnetic moment vectors, one of the following values:
%   * `'auto'`      all moments get the color of the magnetic atom,
%   * `'colorname'` all moments will have the same color defined by the
%                   string, e.g. `'red'`,
%   * `[R G B]`     all moments will have the same color defined by the RGB
%                   code.
% 
% `'scale'`
% : Scaling factor for the lenght of the magnetic moments relative
%   to the length of the shortest bond (if there are no bonds, 3 \\ang 
%   is taken as bond length). Default value is 0.4.
% 
% `'normalize'`
% : If `true`, all moment length will be normalized to the scale
%   factor, default value is `false`.
% 
% `'radius0'`
% : Radius value of arrow body, default value is 0.06 \\ang.
% 
% `'ang'`
% : Angle of the arrow head in degree units, default value is 30 \\deg.
% 
% `'lHead'`
% : Length of the arrow head, default value is 0.5 \\ang.
% 
% `'alpha'`
% : Transparency (alpha value) of the circle, representing the
%   rotation plane of the moments, default value is 0.07.
% 
% #### General paraters
%
% These parameters have the same effect on any of the `swplot.plot...`
% functions.
%
% `'obj'`
% : [spinw] object.
% 
% `'unit'`
% : Unit in which the plotting range is defined. It can be one of the
%   following strings:
%   * `'lu'`        plot range is defined in lattice units (default),
%   * `'xyz'`       plot range is defined in the $xyz$ Cartesian coordinate
%                   system in \\ang units.
%
% `'range'`
% : Defines the plotting range. Depending on the `unit` parameter, the
%   given range can be in lattice units or in units of the $xyz$ Cartesian
%   coordinate system. It is either a matrix with dimensions of $[3\times
%   2]$ where the first and second columns define the lower and upper plot
%   limits respectively. It can be alternatively a vector with three
%   elements `[a,b,c]` which is equivalent to `[0 a;0 b;0 c]`. Default
%   value is `[0 1;0 1;0 1]` to show a single cell.
% 
% `'figure'`
% : Handle of the swplot figure. Default value is the active figure handle.
% 
% `'legend'`
% : Whether to show the plot on the legend, default value is `true`.
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
% `nMesh`
% : Mesh of the ellipse surface, a triangulation class object or an
%   integer that used to generate an icosahedron mesh with $n_{mesh}$
%   number of additional subdivision into triangles. Default value is
%   stored in `swpref.getpref('nmesh')`, see also [swplot.icomesh].
% 
% `nPatch`
% : Number of vertices on any patch object that is not the icosahedron,
%   default value is stored in `swpref.getpref('npatch')`.
%
% ### Output Arguments
% 
% `hFigure`
% : Handle of the swplot figure.
%
% The name of the objects are `'mag'`.
% To find the handles and the corresponding data on these objects, use e.g.
% sObject = swplot.findobj(hFigure,'name','mag')`.
%


% default values
pref = swpref;
fontSize0 = pref.fontsize;
nMesh0    = pref.nmesh;
nPatch0   = pref.npatch;

inpForm.fname  = {'range' 'legend' 'label' 'dtext' 'fontsize' 'radius0' };
inpForm.defval = {[]      true     true    0.1     fontSize0  0.06      };
inpForm.size   = {[-1 -2] [1 1]    [1 1]   [1 1]   [1 1]      [1 1]     };
inpForm.soft   = {true    false    false   false   false      false     };

inpForm.fname  = [inpForm.fname  {'mode' 'color' 'nmesh' 'npatch' 'ang' }];
inpForm.defval = [inpForm.defval {'all'  'auto'  nMesh0  nPatch0  30    }];
inpForm.size   = [inpForm.size   {[1 -4] [1 -5]  [1 1]   [1 1]    [1 1] }];
inpForm.soft   = [inpForm.soft   {false  false   false   false    false }];

inpForm.fname  = [inpForm.fname  {'figure' 'obj' 'unit' 'tooltip' 'copy'}];
inpForm.defval = [inpForm.defval {[]       []    'lu'   true      false }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1] [1 -6] [1 1]     [1 1] }];
inpForm.soft   = [inpForm.soft   {true     true  false  false     false }];

inpForm.fname  = [inpForm.fname  {'shift' 'replace' 'scale' 'normalize' }];
inpForm.defval = [inpForm.defval {[0;0;0] true      0.4      false      }];
inpForm.size   = [inpForm.size   {[3 1]   [1 1]     [1 1]    [1 1]      }];
inpForm.soft   = [inpForm.soft   {false   false     false    false      }];

inpForm.fname  = [inpForm.fname  {'lHead' 'alpha' 'centered' 'translate' 'zoom'}];
inpForm.defval = [inpForm.defval {0.5      0.07   false      false        false}];
inpForm.size   = [inpForm.size   {[1 1]   [1 1]   [1 1]      [1 1]        [1 1]}];
inpForm.soft   = [inpForm.soft   {false   false   false      false        false}];

param = sw_readparam(inpForm, varargin{:});

% find swplot figure
if isempty(param.figure)
    if isempty(param.obj)
        try
            hFigure  = swplot.activefigure;
        catch msg
            warning(msg.message)
            error('plotmag:NoObj','The figure does not contain a SpinW object, use spinw.plot first!')
        end
    else
        hFigure  = swplot.activefigure('plot');
    end
else
    hFigure = param.figure;
end

% takes care of spinw object saved/loaded in/from figure
if isempty(param.obj)
    if isappdata(hFigure,'obj')
        obj = getappdata(hFigure,'obj');
    else
        error('plotmag:NoObj','The figure does not contain a SpinW object, use spinw.plot first!')
    end
else
    if param.copy
        setappdata(hFigure,'obj',copy(param.obj));
    else
        setappdata(hFigure,'obj',param.obj);
    end
    obj = param.obj;
    setappdata(hFigure,'base',obj.basisvector);
end

if obj.symbolic
   warning('plotmag:symbolic','A symbolic magnetic structure can not be plotted.')
   varargout = cell(1:nargout);
   return
end


% the basis vectors in columns.
BV = obj.basisvector;

% set figure title
set(hFigure,'Name', 'SpinW: Magnetic structure');

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
    error('plotmag:WrongInput','The given plotting range is invalid!');
end

switch param.unit
    case 'lu'
        rangelu = [floor(range(:,1)) ceil(range(:,2))];
    case 'xyz'
        % corners of the box
        corners = BV\[range(1,[1 2 2 1 1 2 2 1]);range(2,[1 1 2 2 1 1 2 2]);range(3,[1 1 1 1 2 2 2 2])];
        rangelu = [min(corners,[],2) max(corners,[],2)];
        rangelu = [floor(rangelu(:,1)) ceil(rangelu(:,2))];
    otherwise
        error('plotmag:WrongInput','The given unit string is invalid!');
end

% magnetic supercell for plotting
nExtPlot = diff(rangelu,1,2)'+1;
luOrigin = rangelu(:,1)';
% magnetic moment data
magStr   = obj.magstr('nExt',nExtPlot,'origin',luOrigin,'exact',false);
M        = magStr.S;

% magnetic atom position data
mAtom    = obj.matom;
%nMagAtom = size(mAtom.r,2);

% atom positions in the plotting range
mAtomExt = sw_extendlattice(nExtPlot,mAtom);

pos = bsxfun(@plus,bsxfun(@times,mAtomExt.RRext,nExtPlot'),luOrigin');

% select atom type to plot
switch param.mode
    case 'all'
        % do nothing
    case 'circle'
    case 'arrow'
    case 'none'
    otherwise
        error('plotmag:WrongInput','The given mode string is invalid!');
end

if ~strcmp(param.mode,'none')
    
    % number of unit cells
    nCell = prod(nExtPlot);
    
    % keep track of types of atoms
    aIdx  = repmat(mAtom.idx,[nCell 1])';
    a2Idx = repmat(1:numel(mAtom.idx),[nCell 1])';
    pos  = reshape(pos,3,[]);
    
    % cut out the atoms that are out of range
    switch param.unit
        case 'lu'
            % L>= lower range, L<= upper range
            pIdx = all(bsxfun(@ge,pos,range(:,1)-10*eps) & bsxfun(@le,pos,range(:,2)+10*eps),1);
        case 'xyz'
            % convert to xyz
            posxyz = BV*pos;
            pIdx = all(bsxfun(@ge,posxyz,range(:,1)-10*eps) & bsxfun(@le,posxyz,range(:,2)+10*eps),1);
    end
    
    if ~any(pIdx)
        warning('plotatom:EmptyPlot','There are no magnetic moments in the plotting range!')
        return
    end
    
    M     = M(:,pIdx);
    pos   = pos(:,pIdx);
    aIdx  = aIdx(pIdx);
    a2Idx = a2Idx(pIdx);
    
    % normalization
    if param.normalize
        % normalize moments
        M = bsxfun(@rdivide,M,sqrt(sum(M.^2,1)));
    end
    
    % save magnetic moment vector values into data before rescaling
    MDat = mat2cell([M;pos;a2Idx],7,ones(1,size(M,2)));
    
    % scale moments
    % get the length of the shortest bond
    if ~isempty(obj.coupling.atom2)
        apos1 = obj.matom.r(:,obj.coupling.atom1(1));
        apos2 = obj.matom.r(:,obj.coupling.atom2(1))+double(obj.coupling.dl(:,1));
        lBond = norm(BV*(apos2-apos1));
    else
        lBond = 3;
    end
    
    % normalize the longest moment vector to scale*(shortest bond length)
    if isa(M,'sym')
        warning('sw_plotmag:symcalc','This is a symbolic calculaiton, normalisation can''t be achieved.')
        temp = cellfun(@(x) M./sqrt(sum(x,1).^2),mat2cell(M,3,[1 1 1 1]),'UniformOutput',false);
        [m, ind] = min(cellfun(@(x) sum(x(:)),temp));
        M = temp{ind}*param.scale*lBond;
    else
        M = M/sqrt(max(sum(M.^2,1)))*param.scale*lBond;
    end
    if param.centered
        % double the length for centered moments
        M = 2*M;
    end
    
    % convert to lu units
    Mlu = BV\M;
    
    if param.centered
        % center on atoms
        vpos = cat(3,pos-Mlu/2,pos+Mlu/2);
    else
        vpos = cat(3,pos,pos+Mlu);
    end
    
    % color
    if strcmp(param.color,'auto')
        color = double(obj.unit_cell.color(:,aIdx));
    else
        color = swplot.color(param.color);
    end
    
    % shift positions
    vpos = bsxfunsym(@plus,vpos,BV\param.shift);
    
    % prepare legend labels
    mAtom.name = obj.unit_cell.label(mAtom.idx);
    lLabel = repmat(mAtom.name,[nCell 1]);
    lLabel = lLabel(pIdx);
    
    % plot moment vectors
    swplot.plot('type','arrow','name','mag','position',vpos,'text','',...
        'figure',hFigure,'legend',false,'color',color,'R',param.radius0,...
        'fontsize',param.fontsize,'tooltip',false,'replace',param.replace,...
        'data',MDat,'label',lLabel,'nmesh',param.nmesh,'ang',param.ang,...
        'lHead',param.lHead,'translate',param.translate,'zoom',param.zoom);
    
else
    % don't plot anything just remove previous plot
end


% save range
setappdata(hFigure,'range',struct('range',range,'unit',param.unit));

if nargout > 0
    varargout{1} = hFigure;
end

if param.tooltip
    swplot.tooltip('on',hFigure);
end

end