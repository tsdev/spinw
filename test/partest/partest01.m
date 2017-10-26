function result = partest01(nMat,numWorker)

if nargin == 0
    nMat = 1e3;
end

nMat = round(nMat/numWorker)*numWorker;

% setup
swpref.setpref('usemex',false,'tid',0,'fid',0);

yig = yig_create;
Q = rand(3,nMat);

nSlice  = 4;
nThread = 2;

setenv('OMP_NUM_THREADS',num2str(nThread));

% runs without parallel pool
evalc('delete(gcp(''nocreate''))');
result        = measfun(@spinwavefast_duc,  {yig Q 'optmem' nSlice},false,nMat);
result(end+1) = measfun(@spinwavefast_duc,  {yig Q 'optmem' nSlice},true, nMat);
result(end+1) = measfun(@spinwavefast,      {yig Q 'optmem' nSlice},false,nMat);
result(end+1) = measfun(@spinwavefast,      {yig Q 'optmem' nSlice},true, nMat);
result(end+1) = measfun(@spinwave,          {yig Q 'optmem' nSlice},false,nMat);
result(end+1) = measfun(@spinwave,          {yig Q 'optmem' nSlice},true, nMat);

% run with parpool
evalc(['parpool(' num2str(numWorker) ')']);
result(end+1) = measfun(@spinwavefast,      {yig Q 'optmem' nSlice},false,nMat);
result(end+1) = measfun(@spinwave_spmd,     {yig Q 'optmem' nSlice},false,nMat);

% stop pool
evalc('delete(gcp(''nocreate''))');

end