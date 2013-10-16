function sw_release(verNum, tempDir)
% SW_RELEASE creates a release and packs it into a .zip file.
%
% SW_RELEASE(verNum, {tempDir})
%
% verNum    Version number either  string or a number.
% tempDir   Directory where a temporary subdirectory will be created to
%           store the modified .m files and the final .zip file. By default
%           it is current folder.
% 

if nargin == 0
    help sw_release;
    return
end

% get latest revision number
aDir = pwd;
cd(sw_rootdir);

[statSys, revNum] = system('svn info |grep Revision: |cut -c11-');

if ~statSys
    revNum = str2double(revNum);
else
    revNum = 1;
end

cd(aDir);

if isnumeric(verNum)
    verNum = num2str(verNum);
end

% includes the following comments to every .m file
newLine = sprintf('\n');
revText{1} = ['% $Name: SpinW$ ($Version: ' verNum '$)' newLine];
revText{2} = ['% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)' newLine];
revText{3} = ['% $Revision: ' num2str(revNum) ' $ ($Date: ' date ' $)' newLine];
revText{4} = ['% $License: GNU GENERAL PUBLIC LICENSE$' newLine];
revText{5} = newLine;

% use current directory where the temp directory is created in case no
% tempDir defined
if nargin == 1
    tempDir = pwd;
end

% create temp directory
idx = 1;
tempDirName = [tempDir filesep  'swTemp' num2str(idx,'%02d')];
[statDir, messText] = mkdir(tempDirName);
while ~isempty(messText)
    idx = idx + 1;
    tempDirName = [tempDir filesep  'swTemp' num2str(idx,'%02d')];
    [statDir, messText] = mkdir(tempDirName);
end

if ~statDir
    error('sw_release:CannotCreateDir',['No write access to the ' tempDir ' folder!']);
end

mkdir([tempDirName filesep 'spinw']);
tempDirName0 = tempDirName;
tempDirName = [tempDirName0 filesep 'spinw'];

% copy all files from sw_rootdir to the temp folder
copyfile([sw_rootdir '*'],tempDirName);

% include extra comment to all m files
mFiles = rdir([tempDirName filesep 'm_files' filesep '**' filesep '*.m']);

% go through all files and add comments
for ii = 1:numel(mFiles)
    mLines = {};
    
    fid = fopen(mFiles(ii).name);
    mLines{end+1} = fgets(fid); %#ok<*AGROW>
    mLines{end+1} = fgets(fid);
    
    while numel(mLines{end})>0 && (strcmp(mLines{end}(1),'%'))
        mLines{end+1} = fgets(fid);
    end
    
    % add revision number
    mLines(end+(1:numel(revText))) = revText;
    
    % add remaining lines
    while ~feof(fid)
        mLines{end+1} = fgets(fid);
    end
    fclose(fid);
    
    % open file for rewriting it
    fid = fopen(mFiles(ii).name,'wt');
    
    for jj = 1:numel(mLines)
        fprintf(fid,'%s',mLines{jj});
    end
    fclose(fid);
end

cd(tempDirName0);

zipName = ['spinw_rev' num2str(revNum)];
% compress files
system(['zip -r ' zipName ' * -x \*.DS_Store -x \*.svn']);

zipList = rdir('**/*');

zipDel = {};

for ii = 1:numel(zipList)
    if any(strfind(zipList(ii).name,[filesep '.']))
       zipDel{end+1} = zipList(ii).name;
    end
end

% remove .svn files
system(['zip -d ' zipName ' "*/.*"']);

movefile([zipName '.zip'], aDir);
cd(aDir);
rmdir(tempDirName0,'s');



end
