function [appPath, appName] = sw_apppath()
% returns the path to the app when run in deployed mode

appPath = '';
appName = 'pyspinw';

switch computer
    case 'MACI64'
        %[~, result]      = system(['top -n100 -l1 | grep ' appName ' | awk ''{print $1}''']);
        %[~, result]      = system(['ps -ax | grep "MacOS/' appName ' " | grep -v grep | awk ''{print $1}''']);
        [~, result]      = system(['ps -ax | grep "MacOS/' appName ' " | grep -v grep']);
        exePath = regexp(result,'^.*? (/.*?) ','tokens');
        exePath = exePath{1}{1};
        
        %pid              = min(sscanf(result,'%f'));
        %[status, result] = system(['ps -o comm= -p ' num2str(pid)]);
        %exePath          = strtrim(result);
        
        % get the app path
        appPath = '';
        idx1 = strfind(exePath,[appName '.app']);
        
        if ~isempty(idx1)
            appPath = exePath(1:idx1-2);
        end
        
        appName = [appName '.app'];
    case 'PCWIN64'
        appPath = pwd;
        appName = [appName '.exe'];
    case 'GLNXA64'
        appPath = pwd;
    otherwise
        error('sw_apppath:WrongSystem','Unsupported operating system!')
end

if isempty(appPath)
    error('sw_apppath:MissingAppPath','The app path could not be determined, file an issue on GitHub!')
end


end