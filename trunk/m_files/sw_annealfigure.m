function hFigure = sw_annealfigure
% creates a figure for displaying the status of the annealing simulation
%
% hFigure = SW_ANNEALFIGURE()
%
% See also SW.ANNEAL.
%

% Position of the new figure window.
posFig = get(0,'DefaultFigurePosition');
posFig = [posFig(1:2) 400 400];

% Create new figure.
hFigure = figure;
set(0,'Showhidden','on')

set(hFigure,...
    'Position',      posFig,...
    'Name',          sprintf('Figure %d: SpinW : Annealing status',hFigure),...
    'NumberTitle',   'off',...
    'DockControls',  'off',...
    'PaperType',     'A4',...
    'Tag',           'sw_anneal',...
    'Toolbar',       'figure');

end