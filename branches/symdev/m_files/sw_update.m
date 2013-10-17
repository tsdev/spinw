function sw_update(installDir)
% SW_UPDATE updates the SpinW installation from the internet
%
% SW_UPDATE(installDir)
%
% sw_update creates a new folder with the latest release beside the current
% SpinW installation and add the new version to the working path (and
% removing the old one).
%
% installDir    Folder name, where the new version is installed. Default is
%               the parent folder of the current version of SpinW. If
%               installDir == '.' update will be installed to current
%               folder.
%

% check current version
swVer = sw_version;

% base url, where the sw_download_info file stored
baseUrl = 'https://docs.google.com/uc?export=download&id=0BzFs7CQXhehSRXpjT0dndDNxNUE';


if ~isfield(swVer,'Version')
    answer = getinput('This is a not yet released version of SpinW, update is not recommended! Do you want to continue? (y/n)','yn');
    
    if answer == 'n'
        disp('SpinW update process cancelled!');
        return
    end
    swVer.Revision = 0;
end

if nargin == 0
    installDir = sw_rootdir;
    strIdx = strfind(installDir,filesep);
    installDir = installDir(1:strIdx(end-1));
else
    if installDir(1) == '.'
        installDir = [pwd filesep];
    elseif installDir(end) ~= filesep
        installDir = [installDir filesep];
    end
end


% download the link to the newest version & comments!
% the file format:
%   link to new download
%   revision number
%   message in the next few lines
newInfo = urlread(baseUrl);
% separate lines of text
newInfo = textscan(newInfo, '%s', 'delimiter', sprintf('\n'));
newInfo = newInfo{1};

newLink = newInfo{1};
newRev  = str2double(newInfo{2});
if numel(newInfo)>2
    newMsg  = newInfo(3:end);
else
    newMsg = {};
end

% check whether the online version is newer (compare revision numbers)

fprintf('Current version has a revision number: %d\n',swVer.Revision);
fprintf('New version has a revision number:     %d\n',newRev);

answer = getinput('Do you want to continue? (y/n)','yn');

if answer == 'n'
    disp('SpinW update process cancelled!');
    return
end

fprintf('New version will be installed to: %s\n',installDir);
answer = getinput('Do you want to continue? (y/n)','yn');

if answer(1) == 'n'
    disp('SpinW update process cancelled!');
    return
end

% save new update as a zip file
updateName = 'spinw_update_files.zip';
fprintf('Downloading update from %s... ',newLink);
urlwrite(newLink,[installDir updateName]);
fprintf('ready!\n');

% decompress zip file
zipList = unzip([installDir updateName],installDir);

% get folder name
folName = [installDir strtok(zipList{1}(numel(installDir)+1:end),filesep)];

% remove old SpinW installation from path
disp('Removing path to old SpinW installation!')
rmpath(genpath(sw_rootdir));

% adding new path
fprintf('Adding path to new SpinW installation: %s!\n',folName);
addpath(genpath(folName));

if numel(newMsg)>0
    disp('Release information:')
    disp(repmat('-',[1 60]))
    for ii = 1:numel(newMsg)
        fprintf('\t%s\n',newMsg{ii});
    end
    disp(repmat('-',[1 60]))
end

disp('Removing unnecessary files... ')
delete([installDir updateName]);

disp('In oder to load the new class definitions, issue a ''clear classes'' command before using SpinW!');
disp('In order to reach SpinW after every restart of Matlab, add the following line to your startup.m file:');
fprintf('addpath(genpath(''%s''));\n',folName);
disp('Update was successful!')

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