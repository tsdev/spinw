function sw_compile(varargin)
% compile SpinW code with zeroMQ
%
% SW_COMPILE('option1',value1,...)
%
% Options:
%
% sourcepath        Path to the source files, default value is the output
%                   of sw_rootdir().
% zmqlib            Name of the zeroMQ library file with full path.
%

% To change the icon modify:
% /Applications/MATLAB_R2017a.app/toolbox/compiler/Resources/
%

% location on OSX
switch computer
    case 'MACI64'
        zmqlib0 = '/usr/local/lib/libzmq/libzmq.dylib';
        libname = 'libzmq.dylib';
        % kill all running apps
        [~,~] = system('killall -9 pyspinw');
    case 'PCWIN64'
        zmqlib0 = '';
        libname = 'libzmq.dll';
    case 'GLNXA64'
        zmqlib0 = '';
        libname = 'libzmq.so';
        [~,~] = system('killall -9 pyspinw');
    otherwise
        error('sw_compile:UnsupportedPlatform','Unsupported platform!')
end

swr0    = fileparts(sw_rootdir);

inpForm.fname  = {'sourcepath' 'zmqlib'};
inpForm.defval = {swr0         zmqlib0 };
inpForm.size   = {[1 -1]       [1 -2]  };

param = sw_readparam(inpForm, varargin{:});

fs     = filesep;
swr    = param.sourcepath;

disp('Compiling SpinW...')
tic

if libisloaded('libzmq')
    unloadlibrary('libzmq');
end

% create a temp folder under dev
tPath = [swr fs 'dev' fs 'temp'];
mkdir(tPath)

% copy the zeroMQ library to the temp folder
copyfile(param.zmqlib,[tPath fs libname]);
copyfile([swr fs 'dev' fs 'transplant' fs 'transplantzmq.h'],tPath);
% add a help() function that replaces the Matlab built-in
copyfile([swr fs 'dev' fs 'sw_help.m'],[tPath fs 'help.m']);

% generate the zmqlib header m-file
pwd0 = pwd;
cd(tPath);
% loadlibrary(libname,@libzmq_m,'alias','libzmq')
loadlibrary(libname,'transplantzmq.h','alias','libzmq','mfilename','libzmq_m')
cd(pwd0);

% delete previously compiled files
try %#ok<TRYNC>
    rmdir([sw_rootdir 'dev' fs 'pyspinw' fs 'pyspinw.app'],'s')
    delete([swr fs 'dev' fs 'pyspinw' fs '*'])
end

% add the necessary extra path
addpath([swr fs 'dev'])
addpath([swr fs 'dev' fs 'transplant'])

% compile the code
mccCommand = ['mcc -m '...
    '-d ' swr fs 'dev' fs 'pyspinw '...
    '-o pyspinw '...
    '-a ' swr fs 'dat_files' fs '* '...
    '-a ' swr fs 'external '...
    '-a ' swr fs 'swfiles '...
    '-a ' swr fs 'dev' fs 'waitforgui.m '...
    '-a ' swr fs 'dev' fs 'sw_apppath.m '...
    '-a ' swr fs 'dev' fs 'transplant '...
    '-a ' tPath ' '...
    swr fs 'dev' fs 'pyspinw.m'];
%    '-a ' param.zmqlib ' '...
eval(mccCommand);

% delete the temp directory
rmdir(tPath,'s');

% copy the swfiles into the app
if ismac
    mkdir([swr fs 'dev' fs 'pyspinw' fs 'pyspinw.app' fs 'Source' fs 'swfiles'])
    copyfile([swr fs 'swfiles'],[swr fs 'dev' fs 'pyspinw' fs 'pyspinw.app' fs 'Source' fs 'swfiles'])
    copyfile([swr fs 'dev' fs 'transplant' fs 'transplant_remote.m'],[swr fs 'dev' fs 'pyspinw' fs 'pyspinw.app' fs 'Source'])
    copyfile([swr fs 'dev' fs 'pyspinw.m'],[swr fs 'dev' fs 'pyspinw' fs 'pyspinw.app' fs 'Source'])
    % add new icon
    copyfile([swr '/dev/spinw3.icns'],[swr '/dev/pyspinw/pyspinw.app/Contents/Resources/membrane.icns'],'f');
    system(['touch ' swr '/dev/pyspinw/pyspinw.app']);
else
    warning('sw_compile:MissingFiles','The source files won''t be stored in the copiled app!')
end

% remove unnecessary files
toDel = {'mccExcludedFiles.log' 'readme.txt' 'requiredMCRProducts.txt' 'run_pyspinw.sh'};
for ii = 1:numel(toDel)
    delete([swr fs 'dev' fs 'pyspinw' fs toDel{ii}])
end

% add the script file to the app
copyfile([swr fs 'dev' fs 'pyspinw.sh'],[swr fs 'dev' fs 'pyspinw' fs 'pyspinw.app'])

disp('Done!')
toc

end