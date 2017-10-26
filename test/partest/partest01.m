function result = partest01(nMat,numWorker,nThread,nRun)
% result = partest01(nMat,numWorker,nThread)
%

if nargin == 0
    nMat = 1e3;
end

nMat = round(nMat/numWorker)*numWorker;

% setup
swpref.setpref('usemex',false,'tid',0,'fid',0);

yig = yig_create;
Q = rand(3,nMat);

nSlice  = 4;

setenv('OMP_NUM_THREADS',num2str(nThread));

% runs without parallel pool
evalc('delete(gcp(''nocreate''))');
result        = measfun(@spinwavefast_duc,  {yig Q},false,nSlice,nRun);
result(end+1) = measfun(@spinwavefast_duc,  {yig Q},true, nSlice,nRun);
result(end+1) = measfun(@spinwavefast,      {yig Q},false,nSlice,nRun);
result(end+1) = measfun(@spinwavefast,      {yig Q},true, nSlice,nRun);
result(end+1) = measfun(@spinwave,          {yig Q},false,nSlice,nRun);
result(end+1) = measfun(@spinwave,          {yig Q},true, nSlice,nRun);

% run with parpool
evalc(['parpool(' num2str(numWorker) ')']);
result(end+1) = measfun(@spinwavefast,          {yig Q},false,nSlice,nRun);
result(end+1) = measfun(@spinwave_spmd,         {yig Q},false,nSlice,nRun);
result(end+1) = measfun(@spinwavefast_duc_spmd, {yig Q},false,nSlice,nRun);

% stop pool
evalc('delete(gcp(''nocreate''))');

end