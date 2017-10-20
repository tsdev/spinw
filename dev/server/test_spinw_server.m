%% parameters

folder  = '~/Documents/temp/cache';
nWorker = 0;
portNum = 13001;
jobID   = 'Wa4g3Gj32Z4';

%% produce data

argin   = {sw_model('triAF',1) {[0 0 0] [1 1 0] 1e5} 'fid' 0 'tid' 0};
fun     = 'spinwave';
nargout = 1;

% save .mat file of model
save([folder filesep 'in_' jobID '.mat'],'argin','fun','nargout');

%% run the server
clc
spinw_server(folder,nWorker,portNum);

%% read result

result = load([folder filesep 'in_' jobID '.mat']);

%% direct server

t = tcpip('0.0.0.0',portNum,'NetworkRole','server');
fopen(t);


%%
while 1
    if t.BytesAvailable>0
        char(fread(t,t.BytesAvailable))
    end
end



%% parameters

folder  = '~/temp/cache';
nWorker = 0;
portNum = 13001;
jobID   = 'Wa4g3Gj32Z4';

%% produce data

% empty directory
try
    rmdir(folder,'s')
end
mkdir(folder);
jobIDv = {'J001' 'J002' 'J003'};

% JOB1
% create jobs
argin   = {sw_model('triAF',1) {[0 0 0] [1 1 0] 1e5} 'fid' 1 'tid' 0};
fun     = 'spinwave';
nargout = 1;
% save .mat file of model
save([folder filesep 'in_' jobIDv{1} '.mat'],'argin','fun','nargout');

% JOB2
argin   = {sw_model('squareAF',1) {[0 0 0] [1 0 0] [1 1 0] [0 0 0] 1e5} 'fid' 1 'tid' 0};
fun     = 'spinwave';
nargout = 1;
% save .mat file of model
save([folder filesep 'in_' jobIDv{2} '.mat'],'argin','fun','nargout');

% JOB3
argin   = {sw_model('triAF',1) linspace(0,4,101) 'Evect' linspace(0,8,501) 'nRand' 1e5 'fid' 1 'tid' 0};
fun     = 'powspec';
nargout = 1;
% save .mat file of model
save([folder filesep 'in_' jobIDv{3} '.mat'],'argin','fun','nargout');

% start the server
go swserver
clipboard('copy',['./spinw_server.app/spinw_server.sh ' folder ' 0 ' num2str(portNum)]);
go(folder)

%% start the client

swServer = tcpip('0.0.0.0',portNum, 'NetworkRole', 'client');
fopen(swServer);

fprintf(swServer,['EXEC ' jobIDv{1} ' 1.23:'])
fprintf(swServer,['EXEC ' jobIDv{2} ' 1.23:'])
fprintf(swServer,['EXEC ' jobIDv{3} ' 1.23:'])
pause(1)
fprintf(swServer,['STOP ' jobIDv{3} ':'])


