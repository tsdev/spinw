%% parameters

localFolder  = '~/Documents/temp/cache';
remoteFolder = '/home/lvi05884/spinw_server';
server       = 'lvi05884@isiscompute.nd.rl.ac.uk';
%nWorker      = 32;
portNum      = 13001;

%% produce data

% empty directory
try
    rmdir(localFolder,'s')
end
mkdir(localFolder);
jobIDv = {'J001' 'J002' 'J003' 'J004'};

% JOB1
% create jobs
argin   = {sw_model('triAF',1) {[0 0 0] [1 1 0] 1e5} 'fid' 1 'tid' 0};
fun     = 'spinwave';
nargout = 1;
% save .mat file of model
save([localFolder filesep 'in_' jobIDv{1} '.mat'],'argin','fun','nargout');

% JOB2
argin   = {sw_model('squareAF',1) {[0 0 0] [1 0 0] [1 1 0] [0 0 0] 1e5} 'fid' 1 'tid' 0};
fun     = 'spinwave';
nargout = 1;
% save .mat file of model
save([localFolder filesep 'in_' jobIDv{2} '.mat'],'argin','fun','nargout');

% JOB3
argin   = {sw_model('triAF',1) linspace(0,4,101) 'Evect' linspace(0,8,501) 'nRand' 1e5 'fid' 1 'tid' 0};
fun     = 'powspecfast';
nargout = 1;
% save .mat file of model
save([localFolder filesep 'in_' jobIDv{3} '.mat'],'argin','fun','nargout');

% JOB4
argin   = {yig {Q_N Q_G Q_H 1e4} 'fid' 1 'tid' 0};
fun     = 'spinwavefast';
nargout = 1;
prof    = 1;
% save .mat file of model
save([localFolder filesep 'in_' jobIDv{4} '.mat'],'argin','fun','nargout','prof');

%%
% empty the remote cache
remoteCmd = ['rm -rf ' remoteFolder filesep 'cache' filesep '*'];
system(['ssh ' server ' ''' remoteCmd ''''],'-echo');

% copy the data to the server
system(['rsync -avz ' localFolder filesep '*.mat ' server ':' remoteFolder filesep 'cache'],'-echo');

% forward port
%system(['ssh -L ' num2str(portNum) ':127.0.0.1:' num2str(portNum) ' ' server],'-echo');

%% start remote server + port forwarding

% number of remote threads
[~,nCPU] = system('ssh isis ''./ncpu.sh''');
[~,~,~,~,nCPU] = regexp(nCPU,'\n([0-9]+)\n');
nCore = str2double(nCPU{1}{1})/2;

% set number of workers to the number of cores
nWorker = nCore;

% change the number of workers
go dev
cd server
newCommand = regexprep(fileread('run_server.sh'),' [0-9]+ ',[' ' num2str(nWorker) ' ']);
fid = fopen('run_server.sh','w');
fprintf(fid,newCommand);
fclose(fid);

% start a new terminal window with the server script running
system('open -a Terminal.app run_server.sh')

%% start the client

swServer = tcpip('localhost',portNum, 'NetworkRole', 'client');
fopen(swServer);

fprintf(swServer,['EXEC ' jobIDv{1} ' 1.23:'])
fprintf(swServer,['EXEC ' jobIDv{2} ' 1.23:'])
fprintf(swServer,['EXEC ' jobIDv{3} ' 1.23:'])
fprintf(swServer,['EXEC ' jobIDv{4} ' 1.23:'])
pause(1)
fprintf(swServer,['STOP ' jobIDv{3} ':'])

%% load data from server
system(['rsync -avz ' server ':' remoteFolder filesep 'cache' filesep '*.mat ' localFolder filesep],'-echo');

ii = {};
for ii = 1:numel(jobIDv)
    try
        res{ii} = load([localFolder filesep 'out_' jobIDv{ii} '.mat']);
    end
end