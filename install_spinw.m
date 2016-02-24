function install_spinw()
% installs SpinW
%
% INSTALL_SPINW()
%
% It adds the SpinW folder to the search path and modifies startup.m and
% clears the class definitions. There is a y/n question before every
% operation.
%

% if nargin == 0
%     help sw_install
%     return
% end

% remove old SpinW installation from path
fprintf('\nRemoving path to old SpinW installation if exists!\n')
try %#ok<TRYNC>
    rmpath(genpath(sw_rootdir));
end

% find SpinW folder
folName = fileparts(mfilename('fullpath'));

% adding new path
fprintf('Adding path to new SpinW installation: %s!\n',folName);
ww = warning;
warning('off'); %#ok<WNOFF>
addpath(genpath(folName));
warning(ww);

% remove files aren't needed for new Matlab versions
% functions introduced in R2014a
if ~verLessThan('matlab', '8.1')
    % strjoin()
    fList = dir([folName filesep 'external' filesep 'strjoin*']);
    for ii = 1:numel(fList)
        delete([folName filesep 'external' filesep fList(ii).name]);
    end
    % strsplit
    fList = dir([folName filesep 'external' filesep 'strsplit*']);
    for ii = 1:numel(fList)
        delete([folName filesep 'external' filesep fList(ii).name]);
    end
else
    % rename the functions to be used in Matlab versions prior to R2014a
    % strjoin()
    movefile([folName filesep 'external' filesep 'strjoin0.m'],...
        [folName filesep 'external' filesep 'strjoin.m']);
    % strsplit()
    movefile([folName filesep 'external' filesep 'strsplit0.m'],...
        [folName filesep 'external' filesep 'strsplit.m']);
end

% functions introduced in R2015a
% if ~verLessThan('matlab', '8.5')
%     % uniquetol()
%     fList = dir([folName filesep 'external' filesep 'uniquetol*']);
%     for ii = 1:numel(fList)
%         delete([folName filesep 'external' filesep fList(ii).name]);
%     end
% end

fprintf(['In order to reach SpinW after restarting Matlab, the following\n'...
    'line has to be added to your startup.m file:\n']);
fprintf('  addpath(genpath(''%s''));\n',folName);

% location of Matlab startup file
sfLoc = which('startup');
uPath = userpath;
% remove ':' and ';' characters from the userpath
uPath = [uPath(~ismember(uPath,':;')) filesep 'startup.m'];

% create new startup.m file
if isempty(sfLoc)
    answer = getinput(sprintf(['You don''t have a Matlab startup.m file,\n'...
        'do you want it to be created at %s? (y/n)'],uPath),'yn');
    if answer == 'y'
        fclose(fopen(uPath,'w'));
        sfLoc = uPath;
    end
end

if ~isempty(sfLoc)
    
    answer = getinput(sprintf(['Would you like to add the following line:\n'...
        sprintf('addpath(genpath(''%s''));',folName) '\nto the end of '...
        'your Matlab startup file (%s)? (y/n)'],sfLoc),'yn');
    
    if answer == 'y'
        fid = fopen(sfLoc,'a');
        fprintf(fid,['\n%%###SW_UPDATE\n%% Path to the SpinW installation\n'...
            'addpath(genpath(''%s''));\n%%###SW_UPDATE\n'],folName);
        fclose(fid);
    end
end

answer = getinput(...
    ['\nIn order to refresh the internal class definitions of Matlab (to\n'...
    'access the new SpinW version), issuing the "clear classes" command\n'...
    'is necessary. However this command will also clear all your variables\n'...
    'in the Matlab internal memory. Would you like the updater to issue\n'...
    'the command now, otherwise you can do it manually later.\n'...
    'Do you want to issue the command "clear classes" now? (y/n)'],'yn');

if answer == 'y'
    clear('classes'); %#ok<CLCLS>
    disp('Matlab class memory is refreshed!')
end


disp('The installation of SpinW was successful!')

end

function answer = getinput(message,good)
% get the necessary letter input

answer = ' ';
while ~ismember(answer(1),good)
    answer = input(message,'s');
    if isempty(answer)
        answer = 0;
    end
end
answer = answer(1);

end