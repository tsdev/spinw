function result = measfun(fun,argin,usemex,nMat)
% measure function execution time
%
% result = measfun(fun,argin,usemex,nMat)
%

% check parallel pool
pPool = gcp('nocreate');
if isempty(pPool)
    numWorkers = 1;
else
    numWorkers = pPool.NumWorkers;
end

mexstr = {'nomex' 'mex'};
swpref.setpref('usemex',usemex,'fid',0,'tid',0);
fprintf('Calling %s(), nMat=%d, nWorker=%d, %s...\n',func2str(fun),nMat,numWorkers,mexstr{usemex+1});
tStart = tic;
spec   = fun(argin{:});
tMeas  = toc(tStart);
fprintf('Elapsed time %5.3f s.\n',tMeas);

result.time     = tMeas;
result.fun      = fun;
result.argin    = argin;
result.usemex   = usemex;
result.nmat     = nMat;
result.spec     = spec;

end