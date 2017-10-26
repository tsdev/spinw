function result = measfun(fun,argin,usemex,nMemSlice)
% measure function execution time
%
% result = measfun(fun,argin,usemex,nMat)
%

[~,nThread] = system('echo $OMP_NUM_THREADS');
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

mexstr = {'nomex' 'mex'};
swpref.setpref('usemex',usemex,'fid',0,'tid',0);
fprintf('Calling %s()...\n',func2str(fun));
tStart = tic;
spec   = fun(argin{:});
tMeas  = toc(tStart);
nMemSlice = spec.param.nSlice;
fprintf('... elapsed time %5.3f s, nQ=%d, nWorker=%d, nThread=%s, nMemSlice=%d, %s.\n',...
    tMeas,nQ,nWorker,nThreadStr,nMemSlice,mexstr{usemex+1});

result.time     = tMeas;
result.fun      = fun;
result.argin    = argin;
result.argout   = spec;
result.usemex   = usemex;
result.nQ       = nQ;
result.nWorker  = nWorker;
result.nThread  = nThread;
result.nSlice   = nMemSlice;

end