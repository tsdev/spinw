function result = partest01(nMat,numWorker)

if nargin == 0
    nMat = 1e3;
end

nMat = round(nMat/numWorker)*numWorker;

% setup
swpref.setpref('usemex',false,'tid',0,'fid',0);

yig = yig_create;
Q = rand(3,nMat);

% check parallel pool
pPool = gcp('nocreate');
if isempty(pPool)
    numWorker0 = 1;
else
    numWorker0 = pPool.NumWorkers;
end

if numWorker0~=numWorker
    delete(gcp('nocreate'));
    evalc(['parpool(' num2str(numWorker) ')']);
end


usemex = false;
result = measfun(@spinwavefast,{yig Q},usemex, nMat);
usemex = false;
result(end+1) = measfun(@spinwave_spmd,{yig Q},usemex, nMat);

evalc('delete(gcp)');
% runs without parallel pool
usemex = true;
result(end+1) = measfun(@spinwavefast,{yig Q},usemex, nMat);
usemex = true;
result(end+1) = measfun(@spinwave,{yig Q},usemex, nMat);
usemex = false;
result(end+1) = measfun(@spinwave,{yig Q},usemex, nMat);
usemex = false;
result(end+1) = measfun(@spinwavefast_duc,{yig Q},usemex, nMat);
usemex = true;
result(end+1) = measfun(@spinwavefast_duc,{yig Q},usemex, nMat);

end