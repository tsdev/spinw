%% parameters

folder  = '~/Documents/temp/cache';
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
argin   = {sw_model('triAF',1) {[0 0 0] [1 1 0] 1e3} 'fid' 1 'tid' 0};
fun     = 'spinwave';
nargout = 1;
prof    = 2;

% save .mat file of model
save([folder filesep 'in_' jobIDv{1} '.mat'],'argin','fun','nargout','prof');

% JOB2
argin   = {sw_model('squareAF',1) {[0 0 0] [1 0 0] [1 1 0] [0 0 0] 1e3} 'fid' 1 'tid' 0};
fun     = 'spinwave';
nargout = 1;
prof    = 2;
% save .mat file of model
save([folder filesep 'in_' jobIDv{2} '.mat'],'argin','fun','nargout','prof');

% JOB3
argin   = {sw_model('triAF',1) linspace(0,4,51) 'Evect' linspace(0,8,501) 'nRand' 1e2 'fid' 1 'tid' 0};
fun     = 'powspec';
nargout = 1;
prof    = 2;
% save .mat file of model
save([folder filesep 'in_' jobIDv{3} '.mat'],'argin','fun','nargout','prof');

% start the server
go swserver
spinw_server(folder,0,portNum);

%% start the client

go(folder)
for ii = 1:numel(jobIDv)
    try
        res{ii} = load([folder filesep 'out_' jobIDv{1} '.mat']);
    end
end

