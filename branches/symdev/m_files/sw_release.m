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

swVer = sw_version;

if isfield(swVer,'Version')
    disp('The current version of SpinW is already released!');
    return;
end

% get latest revision number
aDir = pwd;
cd(sw_rootdir);

%[statSys, revNum] = system('svn info |grep Revision: |cut -c11-');
[statSys, revNum] = system('svn up');

revNum = strtrim(revNum);
strIdx = strfind(revNum,' ');
revNum = revNum(strIdx(end):(end-1));

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

% change dots in version name to minus
verNum2 = verNum;
verNum2(verNum2=='.') = '-';

swDirName = ['spinw_' verNum2 '_rev' num2str(revNum)];

mkdir([tempDirName filesep swDirName]);
tempDirName0 = tempDirName;
tempDirName = [tempDirName0 filesep swDirName];

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
    tLines = strtrim(mLines{end});
    
    while numel(tLines)>0 && (strcmp(tLines(1),'%'))
        mLines{end+1} = fgets(fid);
        tLines = strtrim(mLines{end});
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

zipName = [swDirName '.zip'];

% compress files except files and folders with name starting '.' or backup
% files with '~' and 'sw_release.m' file
fList = rdir('**/*');
fListZip = {};

for ii = 1:numel(fList)
    if (~any(strfind(fList(ii).name,[filesep '.']))) && (~any(strfind(fList(ii).name,'~'))) && (~any(strfind(fList(ii).name,'sw_release.m')))
        fListZip{end+1} = fList(ii).name;
    end
end

zip(zipName,fListZip);

movefile(zipName, aDir);
cd(aDir);
rmdir(tempDirName0,'s');

end