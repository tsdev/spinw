function partest(varargin)
% partest01(nQ,numWorker,nThread,nRun)
%

inpForm.fname  = {'nQ'  'nWorker' 'nThread' 'nRun' 'fName'    'nSlice' 'hermit'};
inpForm.defval = {1e2   2         -1        1      'test.mat' 1        false   };
inpForm.size   = {[1 1] [1 -2]    [1 1]     [1 1]  [1 -1]     [1 1]    [1 1]   };

param = sw_readparam(inpForm, varargin{:});

nQ0     = param.nQ;
nWorker = param.nWorker;
nThread = param.nThread;
nRun    = param.nRun;
fName   = param.fName;
nSlice  = param.nSlice;
hermit  = param.hermit;

if nargin == 0
    nQ0 = 1e3;
end

yig = yig_create;
Q = rand(3,nQ0);

% print header
measfun;

% runs without parallel pool
measfun(@spinwavefast_duc,  {yig Q 'hermit', hermit},false,nSlice,nRun,nThread,0,fName);
measfun(@spinwavefast_duc,  {yig Q 'hermit', hermit},true, nSlice,nRun,nThread,0,fName);
measfun(@spinwavefast,      {yig Q 'hermit', hermit},false,nSlice,nRun,nThread,0,fName);
measfun(@spinwavefast,      {yig Q 'hermit', hermit},true, nSlice,nRun,nThread,0,fName);
measfun(@spinwave,          {yig Q 'hermit', hermit},false,nSlice,nRun,nThread,0,fName);
measfun(@spinwave,          {yig Q 'hermit', hermit},true, nSlice,nRun,nThread,0,fName);

for ii = 1:numel(nWorker)
    nQ = round(nQ0/nWorker(ii))*nWorker(ii);
    Q  = rand(3,nQ);

    % run with parpool
    measfun(@spinwavefast,          {yig Q 'hermit', hermit},false,nSlice,nRun,nThread,nWorker(ii),fName);
    measfun(@spinwave_spmd,         {yig Q 'hermit', hermit},false,nSlice,nRun,nThread,nWorker(ii),fName);
    measfun(@spinwavefast_duc_spmd, {yig Q 'hermit', hermit},false,nSlice,nRun,nThread,nWorker(ii),fName);
end
% stop pool
evalc('delete(gcp(''nocreate''))');

end