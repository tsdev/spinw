function hFigure = sw_annealfigure
% hFigure = SW_ANNEALFIGURE() creates a figure for displaying the status of
% the annealing calculation.
%
% See also SW.ANNEAL.
%

if nargin == 0
    help sw_annealfigure;
    return
end

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