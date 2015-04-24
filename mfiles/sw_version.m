function varargout = sw_version()
% returns the installed version of SpinW
%
% SW_VERSION()
%

% read file header from sw.m file
fid = fopen('sw.m');

% first line
fgets(fid);
% skip comments
fLine = strtrim(fgets(fid));
while numel(fLine)>0 && strcmp(fLine(1),'%')
    fLine = strtrim(fgets(fid));
end
% skip empty lines
while numel(fLine)==0 || ~strcmp(fLine(1),'%')
    fLine = strtrim(fgets(fid));
end
% read version information
verLine = {fLine};
while numel(verLine{end})>0 && strcmp(verLine{end}(1),'%')
    verLine{end+1} = strtrim(fgets(fid)); %#ok<*AGROW>
end
verLine = verLine(1:end-1);
fclose(fid);

% read strings enclosed with $ signs
partStr = {};
for ii = 1:numel(verLine)
    verSel = verLine{ii};
    while sum(verSel=='$') > 1
        [~, verSel] = strtok(verSel,'$'); %#ok<*STTOK>
        [partStr{end+1}, verSel] = strtok(verSel,'$');
    end
    
end

nField = numel(partStr);
fieldName = cell(1,nField);
fieldVal = cell(1,nField);

verStruct = struct;

% extract values and save them into a structure
for ii = 1:nField
    [fieldName{ii}, fieldVal{ii}] = strtok(partStr{ii},':');
    fieldVal{ii} = strtrim(fieldVal{ii}(2:end));
    verStruct.(fieldName{ii}) = fieldVal{ii};
end


if nField == 0
    aDir = pwd;
    cd(sw_rootdir);
    [~, revNum] = system('svn info |grep Revision: |cut -c11-');
    cd(aDir);
    revNum = str2double(revNum);
end

% Matlab version & Symbolic Toolbox
v0 = ver;
nSym = strcmp('Symbolic Math Toolbox', {v0.Name});
nSym = find(nSym,1);
if isempty(nSym)
    strSym = 'no Symbolic Math Toolbox installed';
else
    strSym = [v0(nSym).Name ' installed'];
end


if nargout == 0
    if nField == 0
        
        if any(revNum)
            fprintf('This version of SpinW (rev. num. %d) is not released yet!\n',revNum);
        else
            fprintf('This version of SpinW is not released yet!\n');
        end
    else
        disp([verStruct.Name verStruct.Version ' (rev ' num2str(verStruct.Revision) ')']);
        onlineRev = sw_update;
        if onlineRev > str2num(verStruct.Revision) %#ok<ST2NM>
            disp(['Newer version of SpinW is available online (rev. num. ' num2str(onlineRev) '), use the sw_update() function to download it!']);
        else
            disp('You have the latest version of SpinW!')
        end
    end
    fprintf(['MATLAB version: ' version ', ' strSym '\n']);
    
else
    if nField == 0
        if any(revNum)
            varargout{1}.Revision = num2str(revNum);
        else
            varargout{1} = struct;
        end
    else
        varargout{1} = verStruct;
    end
end

end