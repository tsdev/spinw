function compile_macos(varargin)
% compile SpinW code with zeroMQ on MacOS
%
% COMPILE_MACOS('option1',value1,...)
%
% Options:
%
% sourcepath        Path to the source files, default value is the output
%                   of sw_rootdir().
% zmqlib            Name of the zeroMQ library file with full path.
%

if ~ismac
    error('This function works only on MacOS!')
end

% zero-MQ location on OSX
zmqlib0 = '/usr/local/lib/libzmq/libzmq.dylib';
libname = 'libzmq.dylib';
% kill all running apps
[~,~] = system('killall -9 pyspinw');

swr0    = fileparts(sw_rootdir);

inpForm.fname  = {'sourcepath' 'zmqlib'};
inpForm.defval = {swr0         zmqlib0 };
inpForm.size   = {[1 -1]       [1 -2]  };

param = sw_readparam(inpForm, varargin{:});

swr    = param.sourcepath;

disp('Compiling SpinW on MacOS...')
tic

if libisloaded('libzmq')
    unloadlibrary('libzmq');
end

% create a temp folder under dev
cRoot = [swr '/dev/standalone'];
tPath = [cRoot '/MacOS/temp'];
% delete previous temp directory if exists
try %#ok<TRYNC>
    rmdir(tPath,'s');
end
% delete previously compiled files
try %#ok<TRYNC>
    rmdir([cRoot '/MacOS/pyspinw.app'],'s')
    delete([cRoot '/MacOS/*'])
end
mkdir(tPath)

% copy the zeroMQ library to the temp folder
copyfile(param.zmqlib,[tPath filesep libname]);
copyfile([cRoot '/transplant/transplantzmq.h'],tPath);
% add a help() function that replaces the Matlab built-in
copyfile([cRoot '/sw_help.m'],[tPath '/help.m']);

% generate the zmqlib header m-file
pwd0 = pwd;
cd(tPath);
% loadlibrary(libname,@libzmq_m,'alias','libzmq')
loadlibrary(libname,'transplantzmq.h','alias','libzmq','mfilename','libzmq_m')
cd(pwd0);

% add the necessary extra path
addpath([swr '/dev/standalone'])
addpath([cRoot '/transplant'])

% compile the code
mccCommand = ['mcc -m '...
    '-d ' cRoot '/MacOS '...
    '-o pyspinw '...
    '-a ' swr '/dat_files/* '...
    '-a ' swr '/external '...
    '-a ' swr '/swfiles '...
    '-a ' cRoot '/waitforgui.m '...
    '-a ' cRoot '/sw_apppath.m '...
    '-a ' cRoot '/transplant '...
    '-a ' tPath ' '...
    cRoot '/pyspinw.m'];

eval(mccCommand);

% delete the temp directory
rmdir(tPath,'s');

% copy the swfiles into the app
mkdir([cRoot '/MacOS/pyspinw.app/Source/swfiles'])
copyfile([swr '/swfiles'],[cRoot '/MacOS/pyspinw.app/Source/swfiles'])
copyfile([cRoot '/transplant/transplant_remote.m'],[cRoot '/MacOS/pyspinw.app/Source'])
copyfile([cRoot '/pyspinw.m'],[cRoot '/MacOS/pyspinw.app/Source'])
% add new icon
copyfile([cRoot '/icon/spinw3.icns'],[cRoot '/MacOS/pyspinw.app/Contents/Resources/membrane.icns'],'f');
system(['touch ' cRoot '/MacOS/pyspinw.app']);

% remove unnecessary files
toDel = {'mccExcludedFiles.log' 'readme.txt' 'requiredMCRProducts.txt' 'run_pyspinw.sh'};
for ii = 1:numel(toDel)
    delete([cRoot '/MacOS/' toDel{ii}])
end

% add the script file to the app
copyfile([cRoot '/pyspinw_macos.sh'],[cRoot '/MacOS/pyspinw.app/pyspinw.sh'])

disp('Done!')
toc

end