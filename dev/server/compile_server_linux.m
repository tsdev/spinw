function compile_server_linux(varargin)
% compile SpinW server on Linux
%
% COMPILE_SERVER_LINUX('option1',value1,...)
%
% Options:
%
% sourcepath        Path to the source files, default value is the output
%                   of sw_rootdir().
% zmqlib            Name of the zeroMQ library file with full path.
%

%if ~isunix || ismac
%    error('This function works only on Linux!')
%end

% kill all running apps
[~,~] = system('killall -9 spinw_server');

swr0    = fileparts(sw_rootdir);

inpForm.fname  = {'sourcepath'};
inpForm.defval = {swr0        };
inpForm.size   = {[1 -1]      };

param = sw_readparam(inpForm, varargin{:});

swr    = param.sourcepath;

disp('Compiling SpinW Server on Linux...')
tic

% create a temp folder under dev
cRoot = [swr '/dev/server'];
tPath = [cRoot '/Linux/temp'];
try %#ok<TRYNC>
    rmdir(tPath,'s');
end
% delete previously compiled files
try %#ok<TRYNC>
    rmdir([cRoot '/Linux/Source'],'s')
    delete([cRoot '/Linux/*'])
end
mkdir(tPath)

% add the necessary extra path
addpath([swr '/dev/server'])

% compile the code
mccCommand = ['mcc -m '...
    '-d ' cRoot '/Linux '...
    '-o spinw_server '...
    '-a ' swr '/dat_files/* '...
    '-a ' swr '/external '...
    '-a ' swr '/swfiles '...
    '-a ' tPath ' '...
    cRoot '/spinw_server.m'];

eval(mccCommand);

% delete the temp directory
rmdir(tPath,'s');

% remove unnecessary files
toDel = {'mccExcludedFiles.log' 'readme.txt' 'requiredMCRProducts.txt' 'run_spinw_server.sh'};
for ii = 1:numel(toDel)
    delete([cRoot '/Linux/' toDel{ii}])
end

% add the script file to the app
copyfile([cRoot '/spinw_server_linux.sh'],[cRoot '/Linux/spinw_server.sh'])

disp('Done!')
toc

end