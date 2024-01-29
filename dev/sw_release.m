function sw_release(verNum, tempDir)
% creates a SpinW release and packs it into a .zip file
%
% SW_RELEASE()
%
% SW_RELEASE(verNum, {tempDir})
%
% verNum    Version number either  string or a number.
% tempDir   Directory where a temporary subdirectory will be created to
%           store the modified .m files and the final .zip file. By default
%           it is current folder.
%
% See also SW_INITIALIZE, SW_VERSION, SW_ROOTDIR.
%

if nargin == 0
    swhelp sw_release
    return
end

% get version information for the release
swVer = sw_version;

if ~isempty(swVer.Version)
    disp('The current version of SpinW is already released!');
    return
end

if isnumeric(verNum)
    verNum = num2str(verNum);
end

% includes the following comments to every .m file
newLine    = char(10); %#ok<CHARTEN>
revText{1} = ['% $Name: ' swVer.Name '$ ($Version: ' verNum '$)' newLine];
revText{2} = ['% $Author: ' swVer.Author '$ ($Contact: ' swVer.Contact '$)' newLine];
revText{3} = ['% $Revision: ' swVer.Release '$ ($Date: ' date '$)' newLine];
revText{4} = ['% $License: ' swVer.License '$' newLine];
revText{5} = newLine;

% current directory
aDir = pwd;

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

swDirName = ['spinw' verNum2(1) '_R' num2str(swVer.Release)];

mkdir([tempDirName filesep swDirName]);
tempDirName0 = tempDirName;
tempDirName = [tempDirName0 filesep swDirName];

% initialize the symmetry.dat file
sw_initialize;

% copy all files from sw_rootdir to the temp folder
copyfile([sw_rootdir '*'],tempDirName);

% include extra comment to all m files
swFiles = rdir([tempDirName filesep 'swfiles' filesep '**' filesep '*.m']);

% go through all files and add comments
for ii = 1:numel(swFiles)
    mLines = {};
    
    fid = fopen(swFiles(ii).name);
    mLines{end+1} = fgets(fid); %#ok<*AGROW>
    mLines{end+1} = fgets(fid);
    tLines = strtrim(mLines{end});
    
    while numel(tLines)>0 && (strcmp(tLines(1),'%'))
        mLines{end+1} = fgets(fid);
        if ischar(mLines{end})
            tLines = strtrim(mLines{end});
        else
            tLines = 0;
            mLines{end} = newline;
        end
    end
    
    % add release number
    mLines(end+(1:numel(revText))) = revText;
    
    % add remaining lines
    while ~feof(fid)
        mLines{end+1} = fgets(fid);
    end
    fclose(fid);
    
    % open file for rewriting it
    fid = fopen(swFiles(ii).name,'wt');
    
    for jj = 1:numel(mLines)
        fprintf(fid,'%s',mLines{jj});
    end
    fclose(fid);
end

% change the Contents file
cId = fopen([tempDirName filesep 'Contents.m'],'w');
fprintf(cId,['%% SpinW\n%% Version ' verNum ' (R' num2str(swVer.Release) ') ' swVer.Date '\n%%']);
fclose(cId);

cd(tempDirName0);

zipName = [swDirName '.zip'];

% compress files except files and folders with name starting '.' or backup
% files with '~' and 'sw_release.m' file
fList = rdir('**/*');
fListZip = {};
if ispc
    sp = '\';
else
    sp = '/';
end
dirname = [pwd sp];

for ii = 1:numel(fList)
    if (~any(strfind(fList(ii).name,[filesep '.']))) && (~any(strfind(fList(ii).name,'~'))) ...
            && (~any(strfind(fList(ii).name,[filesep 'dev' filesep]))) ...
            && (~any(strfind(fList(ii).name,[filesep 'docs' filesep]))) ...
            && (~any(strfind(fList(ii).name,[filesep 'test' filesep]))) ...
            && (~any(strfind(fList(ii).name,[filesep 'tutorials' filesep]))) ...
            && (~any(strfind(fList(ii).folder,'git')))
        fListZip{end+1} = fList(ii).name;
    end
end

zip(zipName,fListZip);

movefile(zipName, aDir);
cd(aDir);
rmdir(tempDirName0,'s');

end