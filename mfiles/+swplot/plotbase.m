function plotbase(varargin)
% plots the edges of unit cells on swplot figure
%
% SWPLOT.PLOTCELL('Option1',Value1,...)
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
%

% default values
range0 = [0 1;0 1;0 1];
col0   = swplot.color({'red' 'green' 'blue'});
d0     = ones(1,3);

inpForm.fname  = {'range' 'mode'    'figure' 'color' 'R'   'alpha' 'lhead' 'd'   'dtext' 'label'};
inpForm.defval = {range0  'default' []       col0    0.06  30      0.5     d0    0.5     true   };
inpForm.size   = {[-1 -2] [1 -3]    [1 1]    [-4 -5] [1 1] [1 1]   [1 1]   [1 3] [1 1]   [1 1]  };
inpForm.soft   = {false   false     true     false   false false   false   false false   false  };

param = sw_readparam(inpForm, varargin{:});

if isempty(param.figure)
    hFigure  = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

% basis vectors
BV = swplot.base(hFigure);

% convert d from xyz to base
d = param.d/(BV');

R = bsxfun(@minus,cat(3,zeros(3),eye(3)),d');
% plot the arrows
swplot.plot('type','arrow','position',R,'figure',hFigure,'color',param.color,...
    'name','abc','legend',false,'tooltip',false);

if param.label
    % convert d from xyz to base
    Rtext = bsxfun(@minus,bsxfun(@times,eye(3),1+param.dtext./sqrt(sum(BV.^2,1))),d');
    
    swplot.plot('type','text','position',Rtext,'text',{'a' 'b' 'c'},...
        'figure',hFigure,'legend',false,'tooltip',false);
end

end