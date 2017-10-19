function compile_server_macos(varargin)
% compile SpinW code with zeroMQ on MacOS
%
% COMPILE_SERVER_MACOS('option1',value1,...)
%
% Options:
%
% sourcepath        Path to the source files, default value is the output
%                   of sw_rootdir().
%

if ~ismac
    error('This function works only on MacOS!')
end

% kill all running apps
[~,~] = system('killall -9 spinw_server');

swr0    = fileparts(sw_rootdir);

inpForm.fname  = {'sourcepath'};
inpForm.defval = {swr0        };
inpForm.size   = {[1 -1]      };

param = sw_readparam(inpForm, varargin{:});

swr    = param.sourcepath;

disp('Compiling SpinW Server on MacOS...')
tic

% create a temp folder under dev
cRoot = [swr '/dev/server'];
tPath = [cRoot '/MacOS/temp'];
% delete previous temp directory if exists
try %#ok<TRYNC>
    rmdir(tPath,'s');
end
% delete previously compiled files
try %#ok<TRYNC>
    rmdir([cRoot '/MacOS/spinw_server.app'],'s')
    delete([cRoot '/MacOS/*'])
end
mkdir(tPath)

% add the necessary extra path
addpath([swr '/dev/server'])

% compile the code
mccCommand = ['mcc -m '...
    '-d ' cRoot '/MacOS '...
    '-o spinw_server '...
    '-a ' swr '/dat_files/* '...
    '-a ' swr '/external '...
    '-a ' swr '/swfiles '...
    '-a ' tPath ' '...
    cRoot '/spinw_server.m'];

eval(mccCommand);

% delete the temp directory
rmdir(tPath,'s');

% add new icon
copyfile([cRoot '/icon/spinw3.icns'],[cRoot '/MacOS/spinw_server.app/Contents/Resources/membrane.icns'],'f');
system(['touch ' cRoot '/MacOS/spinw_server.app']);

% remove unnecessary files
toDel = {'mccExcludedFiles.log' 'readme.txt' 'requiredMCRProducts.txt' 'run_spinw_server.sh'};
for ii = 1:numel(toDel)
    delete([cRoot '/MacOS/' toDel{ii}])
end

% add the script file to the app
copyfile([cRoot '/spinw_server_macos.sh'],[cRoot '/MacOS/spinw_server.app/spinw_server.sh'])

disp('Done!')
toc

end