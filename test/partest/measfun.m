function result = measfun(fun,argin,usemex,nMemSlice)
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
    
if nargin > 4
    argin = [argin(:) {'optmem' nMemSlice}];
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
tStart = tic;
spec   = fun(argin{:});
tMeas  = toc(tStart);
nMemSlice = spec.param.nSlice;
fprintf('... %6.1f px/sec, elapsed time %5.3f s, nQ=%d, nMode=%d, nWorker=%d, nThread=%s, nMemSlice=%d, %s.\n',...
    nQ/tMeas,tMeas,nQ,nMode,nWorker,nThreadStr,nMemSlice,mexstr{usemex+1});

result.time     = tMeas;
result.fun      = fun;
result.argin    = argin;
result.argout   = spec;
result.usemex   = usemex;
result.nQ       = nQ;
result.nWorker  = nWorker;
result.nThread  = nThread;
result.nSlice   = nMemSlice;
result.pps      = nQ/tMeas;

end