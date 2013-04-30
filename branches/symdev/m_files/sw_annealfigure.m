function hFigure = sw_annealfigure
% hFigure = SW_ANNEALFIGURE() creates a figure for displaying the status of
% the annealing calculation.
%
% See also SW.ANNEAL.
%

% Position of the new figure window.
posFig = [300 100 400 400];

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

uData = [];
set(hFigure,'UserData',uData);

end