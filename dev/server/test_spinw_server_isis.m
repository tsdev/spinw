%% parameters

localFolder  = '~/Documents/temp/cache';
remoteFolder = '/home/lvi05884/spinw_server';
server       = 'lvi05884@isiscompute.nd.rl.ac.uk';
nWorker      = 0;
portNum      = 13001;

%% produce data

% empty directory
try
    rmdir(localFolder,'s')
end
mkdir(localFolder);
jobIDv = {'J001' 'J002' 'J003'};

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
fun     = 'powspec';
nargout = 1;
% save .mat file of model
save([localFolder filesep 'in_' jobIDv{3} '.mat'],'argin','fun','nargout');

% empty the remote cache
remoteCmd = ['rm -rf ' remoteFolder filesep 'cache' filesep '*'];
system(['ssh ' server ' ''' remoteCmd '''&']);

% copy the data to the server
system(['rsync -avz ' localFolder filesep '*.mat ' server ':' remoteFolder filesep 'cache'],'-echo');

% start the server remotely on ISISCOMPUTE
remoteCmd = [remoteFolder filesep 'spinw_server.sh ' remoteFolder filesep 'cache' ' ' num2str(nWorker) ' ' num2str(portNum)];
system(['ssh ' server ' ''' remoteCmd '''&']);

% wait for the server to startup
pause(10);

% forward port
system(['ssh -L ' num2str(portNum) ':127.0.0.1:' num2str(portNum) ' ' server],'-echo');

%% start the client

swServer = tcpip('0.0.0.0',portNum, 'NetworkRole', 'client');
fopen(swServer);

fprintf(swServer,['EXEC ' jobIDv{1} ' 1.23:'])
fprintf(swServer,['EXEC ' jobIDv{2} ' 1.23:'])
fprintf(swServer,['EXEC ' jobIDv{3} ' 1.23:'])
pause(1)
fprintf(swServer,['STOP ' jobIDv{3} ':'])


