function varargout = plotbase(varargin)
% plots the edges of unit cells on swplot figure
%
% SWPLOT.PLOTBASE('Option1',Value1,...)
%
% hFigure = SWPLOT.PLOTBASE('Option1',Value1,...)
%
% Options:
%
% label     Logical variable, plots abc labels if true. Default is true.
% figure    Handle of the swplot figure. Default is the selected figure.
% color     Color of the arrows, default is red-green-blue for abc, stored
%           in the columns of a 3x3 matrix.
% R         Radius value of arrow body, default is 0.06.
% alpha     Head angle of the arrow in degree units, default is 30 degree.
% lHead     Length of the arrow head, default value is 0.5.
% d         Distance from origin in xyz units.
% dtext     Distance from arrow and in xyz units.
% translate If true, all plot objects will be translated to the figure
%           center. Default is true.
% zoom      If true, figure will be automatically zoomed to the ideal size.
%           Default is true.
%

% default values
range0  = [0 1;0 1;0 1];
col0    = swplot.color({'red' 'green' 'blue'});
d0      = ones(1,3);
nMesh0  = swpref.getpref('nmesh',[]);
nPatch0 = swpref.getpref('npatch',[]);

inpForm.fname  = {'range' 'mode'    'figure' 'color' 'R'   'alpha' 'lhead'};
inpForm.defval = {range0  'default' []       col0    0.06  30      0.5    };
inpForm.size   = {[-1 -2] [1 -3]    [1 1]    [-4 -5] [1 1] [1 1]   [1 1]  };
inpForm.soft   = {false   false     true     false   false false   false  };

inpForm.fname  = [inpForm.fname  {'d'   'dtext' 'label' 'tooltip' 'translate'}];
inpForm.defval = [inpForm.defval {d0    0.5     true    true      true       }];
inpForm.size   = [inpForm.size   {[1 3] [1 1]   [1 1]   [1 1]     [1 1]      }];
inpForm.soft   = [inpForm.soft   {false false   false   false     false      }];

inpForm.fname  = [inpForm.fname  {'zoom' 'replace' 'nmesh' 'npatch' 'unit'}];
inpForm.defval = [inpForm.defval {true   true      nMesh0  nPatch0  'lu'  }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]     [1 1]   [1 1]    [1 -6]}];
inpForm.soft   = [inpForm.soft   {false  false     false   false    false }];

param = sw_readparam(inpForm, varargin{:});

if isempty(param.figure)
    hFigure  = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

% range of previous plot
rDat   = getappdata(hFigure,'range');
range0 = rDat.range;
unit   = rDat.unit;

if isempty(range0)
    range0 = [0;0;0];
else
    range0 = range0(:,1);
end

% basis vectors
BV = swplot.base(hFigure);

switch unit
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

R = bsxfun(@minus,cat(3,zeros(3),eye(3)),d);

% generate data
data = repmat({BV},[1 3]);

if param.label
    % convert d from xyz to base
    Rtext = bsxfun(@minus,bsxfun(@times,eye(3),1+param.dtext./sqrt(sum(BV.^2,1))),d);
    
    swplot.plot('type','text','position',Rtext,'text',{'a' 'b' 'c'},'data',data,...
        'figure',hFigure,'legend',false,'tooltip',false,'translate',false,...
        'zoom',false,'name','base_label','replace',param.replace);
end

% plot the arrows
swplot.plot('type','arrow','position',R,'figure',hFigure,'color',param.color,...
    'name','base','legend',false,'tooltip',false,'translate',param.translate,...
    'zoom',param.zoom,'data',data,'replace',param.replace,'npatch',param.npatch);

if param.tooltip
    swplot.tooltip('on',hFigure);
end

if nargout>0
    varargout{1} = hFigure;
end

end