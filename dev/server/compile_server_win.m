function compile_win(varargin)
% compile SpinW code with zeroMQ on Windows
%
% COMPILE_WIN('option1',value1,...)
%
% Options:
%
% sourcepath        Path to the source files, default value is the output
%                   of sw_rootdir().
% zmqlib            Name of the zeroMQ library file with full path.
%

if ~ispc
    error('This function works only on Windows!')
end

% zero-MQ location on Win
zmqlib0 = 'C:\Program Files\ZeroMQ 4.0.4\bin\libzmq-v120-mt-4_0_4.dll';
libname = 'libzmq.dll';
% kill all running apps
[~,~] = system('taskkill.exe /F /IM pyspinw.exe /T');

swr0    = fileparts(sw_rootdir);

inpForm.fname  = {'sourcepath' 'zmqlib'};
inpForm.defval = {swr0         zmqlib0 };
inpForm.size   = {[1 -1]       [1 -2]  };

param = sw_readparam(inpForm, varargin{:});

swr    = param.sourcepath;

disp('Compiling SpinW on Windows...')
tic

if libisloaded('libzmq')
    unloadlibrary('libzmq');
end

% create a temp folder under dev
cRoot = [swr '\dev\standalone'];
tPath = [cRoot '\Win\temp'];
% delete previous temp directory if exists
try %#ok<TRYNC>
    rmdir(tPath,'s');
end
% delete previously compiled files
try %#ok<TRYNC>
    rmdir([cRoot '\Win\pyspinw'],'s')
    delete([cRoot '\Win\*'])
end
mkdir(tPath)

% copy the zeroMQ library to the temp folder
copyfile(param.zmqlib,[tPath filesep libname]);
copyfile([cRoot '\transplant\transplantzmq.h'],tPath);
% add a help() function that replaces the Matlab built-in
copyfile([cRoot '\sw_help.m'],[tPath '\help.m']);

% generate the zmqlib header m-file
pwd0 = pwd;
cd(tPath);
% loadlibrary(libname,@libzmq_m,'alias','libzmq')
loadlibrary(libname,'transplantzmq.h','alias','libzmq','mfilename','libzmq_m')
cd(pwd0);

% add the necessary extra path
addpath([swr '\dev\standalone'])
addpath([cRoot '\transplant'])

% compile the code
mccCommand = ['mcc -m '...
    '-d ' cRoot '\Win '...
    '-o pyspinw '...
    '-a ' swr '\dat_files\* '...
    '-a ' swr '\external '...
    '-a ' swr '\swfiles '...
    '-a ' cRoot '\waitforgui.m '...
    '-a ' cRoot '\sw_apppath.m '...
    '-a ' cRoot '\transplant '...
    '-a ' tPath ' '...
    cRoot '\pyspinw.m'];

eval(mccCommand);

% delete the temp directory
unloadlibrary('libzmq')
rmdir(tPath,'s');

% copy the swfiles into the app
mkdir([cRoot '\Win\Source\swfiles'])
copyfile([swr '\swfiles'],[cRoot '\Win\Source\swfiles'])
copyfile([cRoot '\transplant\transplant_remote.m'],[cRoot '\Win\Source'])
copyfile([cRoot '\pyspinw.m'],[cRoot '\Win\Source'])
% add new icon
%copyfile([cRoot '\icon\spinw3.icns'],[cRoot '\Win\pyspinw.app\Contents\Resources\membrane.icns'],'f');

% remove unnecessary files
toDel = {'mccExcludedFiles.log' 'readme.txt' 'requiredMCRProducts.txt'};
for ii = 1:numel(toDel)
    delete([cRoot '\Win\' toDel{ii}])
end

disp('Done!')
toc

end