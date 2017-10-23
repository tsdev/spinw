%% parameters

localFolder  = '~/Documents/temp/cache';
remoteFolder = '/home/l_mc01/mpi/toth/spinw_server/cache';
nWorker = 0;
portNum = 13001;

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

% copy the data to the server
for ii = 1:numel(jobIDv)
    remoteFile = ['l_mc01@mcc.psi.ch:' remoteFolder filesep 'in_' jobIDv{ii} '.mat'];
    localFile  = [localFolder filesep 'in_' jobIDv{ii} '.mat'];
    system(['rsync -avz ' localFile ' ' remoteFile],'-echo');
end

% start the server remotely on MCC
clipboard('copy',['./spinw_server_linux.sh ' remoteFolder ' 0 ' num2str(portNum)]);


%% start the client

swServer = tcpip('0.0.0.0',portNum, 'NetworkRole', 'client');
fopen(swServer);

fprintf(swServer,['EXEC ' jobIDv{1} ' 1.23:'])
fprintf(swServer,['EXEC ' jobIDv{2} ' 1.23:'])
fprintf(swServer,['EXEC ' jobIDv{3} ' 1.23:'])
pause(1)
fprintf(swServer,['STOP ' jobIDv{3} ':'])


