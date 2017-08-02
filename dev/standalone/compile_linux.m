function compile_linux(varargin)
% compile SpinW code with zeroMQ on Linux
%
% COMPILE_LINUX('option1',value1,...)
%
% Options:
%
% sourcepath        Path to the source files, default value is the output
%                   of sw_rootdir().
% zmqlib            Name of the zeroMQ library file with full path.
%

if ~isunix || ismac
    error('This function works only on Linux!')
end

% library location on Linux
zmqlib0 = '/usr/local/lib64/libzmq.so';
libname = 'libzmq.so';
% kill all running apps
[~,~] = system('killall -9 pyspinw');

swr0    = fileparts(sw_rootdir);

inpForm.fname  = {'sourcepath' 'zmqlib'};
inpForm.defval = {swr0         zmqlib0 };
inpForm.size   = {[1 -1]       [1 -2]  };

param = sw_readparam(inpForm, varargin{:});

swr    = param.sourcepath;

disp('Compiling SpinW on Linux...')
tic

if libisloaded('libzmq')
    unloadlibrary('libzmq');
end

% create a temp folder under dev
cRoot = [swr '/dev/standalone'];
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

% copy the zeroMQ library to the temp folder
system(['cp --preserve=links --remove-destination ' param.zmqlib ' ' tPath filesep libname]);
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
    '-d ' cRoot '/Linux '...
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

mkdir([cRoot '/Linux/Source/swfiles'])
copyfile([swr '/swfiles'],[cRoot '/Linux/Source/swfiles'])
copyfile([cRoot '/transplant/transplant_remote.m'],[cRoot '/Linux/Source'])
copyfile([cRoot '/pyspinw.m'],[cRoot '/Linux/Source'])

% remove unnecessary files
toDel = {'mccExcludedFiles.log' 'readme.txt' 'requiredMCRProducts.txt'};
for ii = 1:numel(toDel)
    delete([cRoot '/Linux/' toDel{ii}])
end

% add the script file to the app
copyfile([cRoot '/pyspinw_linux.sh'],[cRoot '/Linux/pyspinw.sh'])

disp('Done!')
toc

end