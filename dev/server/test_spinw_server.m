%% parameters

folder  = '~/Documents/temp/cache';
nWorker = 0;
portNum = 13001;
jobID   = 'Wa4g3Gj32Z';

%% produce data

argin   = {sw_model('triAF',1) {[0 0 0] [1 1 0] 501}};
fun     = 'spinwave';
nargout = 1;

% save .mat file of model
save([folder filesep 'in_' jobID '.mat'],'argin','fun','nargout');

%% run the server
clc
spinw_server(folder,nWorker,portNum);

%% read result

result = load([folder filesep 'out_' jobID '.mat']);

%% direct server

t = tcpip('0.0.0.0',portNum,'NetworkRole','server');
fopen(t);


%%
while 1
    if t.BytesAvailable>0
        char(fread(t,t.BytesAvailable))
    end
end