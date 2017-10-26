function result = measfun(fun,argin,usemex,nMemSlice,nTest)
% measure function execution time
%
% result = measfun(fun,argin,usemex,nMat)
%

nThread = getenv('OMP_NUM_THREADS');
nThread = str2double(nThread);
if isempty(nThread) || isnan(nThread)
    nThread = -1;
end

% check parallel pool
pPool = gcp('nocreate');
if isempty(pPool)
    nWorker = 1;
else
    nWorker = pPool.NumWorkers;
end

if nThread<1
    nThreadStr = 'max';
else
    nThreadStr = num2str(nThread);
end

if iscell(argin{2})
    argin{2} = sw_qscan(argin{2});
end
nQ = size(argin{2},2);
    
if nargin > 3
    argin = [argin(:)' {'optmem' nMemSlice}];
end

if nargin<5
    nTest = 1;
end

% number of spin wave modes
nMode = argin{1}.nmagext*2;
% number of Q points triples for incommensurate structure
if any(mod(argin{1}.mag_str.k(:,1),1))
    nQ = nQ*3;
end

mexstr = {'nomex' 'mex'};
swpref.setpref('usemex',usemex,'fid',0,'tid',0);
fprintf('Calling %s()...\n',func2str(fun));

for ii = 1:nTest
    tStart    = tic;
    spec      = fun(argin{:});
    tMeas(ii) = toc(tStart);
end

pps   = nQ./tMeas;

if nTest == 1
    tMeas = round(tMeas*1e3)/1e3;
    pps   = round(pps*1e3)/1e3;
end

tMeasStr = err2str(mean(tMeas),std(tMeas));
ppsStr   = err2str(mean(pps),std(pps));

nMemSlice = spec.param.nSlice;
fprintf('... %s px/sec, elapsed time %s s, nQ=%d, nMode=%d, nWorker=%d, nThread=%s, nMemSlice=%d, %s.\n',...
    ppsStr,tMeasStr,nQ,nMode,nWorker,nThreadStr,nMemSlice,mexstr{usemex+1});

result.time     = tMeas;
result.fun      = fun;
result.argin    = argin;
result.argout   = spec;
result.usemex   = usemex;
result.nQ       = nQ;
result.nWorker  = nWorker;
result.nThread  = nThread;
result.nSlice   = nMemSlice;
result.pps      = pps;

end