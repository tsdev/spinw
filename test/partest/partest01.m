function partest01(varargin)
% partest01(nQ,numWorker,nThread,nRun)
%

inpForm.fname  = {'nQ'  'nWorker' 'nThread' 'nRun' 'fName'   };
inpForm.defval = {1e2   2         -1        1      'test.mat'};
inpForm.size   = {[1 1] [1 1]     [1 1]     [1 1]  [1 -1]    };

param = sw_readparam(inpForm, varargin{:});

nQ      = param.nQ;
nWorker = param.nWorker;
nThread = param.nThread;
nRun    = param.nRun;
fName   = param.fName;


if nargin == 0
    nQ = 1e3;
end

nQ = round(nQ/nWorker)*nWorker;

% setup
swpref.setpref('usemex',false,'tid',0,'fid',0);

yig = yig_create;
Q = rand(3,nQ);

nSlice  = 4;

if nThread > 0
    setenv('OMP_NUM_THREADS',num2str(nThread));
else
    setenv('OMP_NUM_THREADS','');
end
% print header
measfun;

% runs without parallel pool
evalc('delete(gcp(''nocreate''))');
measfun(@spinwavefast_duc,  {yig Q},false,nSlice,nRun,fName);
measfun(@spinwavefast_duc,  {yig Q},true, nSlice,nRun,fName);
measfun(@spinwavefast,      {yig Q},false,nSlice,nRun,fName);
measfun(@spinwavefast,      {yig Q},true, nSlice,nRun,fName);
measfun(@spinwave,          {yig Q},false,nSlice,nRun,fName);
measfun(@spinwave,          {yig Q},true, nSlice,nRun,fName);

% run with parpool
evalc(['parpool(' num2str(nWorker) ')']);
measfun(@spinwavefast,          {yig Q},false,nSlice,nRun,fName);
measfun(@spinwave_spmd,         {yig Q},false,nSlice,nRun,fName);
measfun(@spinwavefast_duc_spmd, {yig Q},false,nSlice,nRun,fName);

% stop pool
evalc('delete(gcp(''nocreate''))');

end