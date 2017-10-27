function partest2(varargin)
% partest01(nQ,numWorker,nThread,nRun)
%

inpForm.fname  = {'nQ'  'nWorker' 'nThread' 'nRun' 'fName'    'nSlice' 'hermit'};
inpForm.defval = {1e3   2         -1        1      'test.mat' 1        false   };
inpForm.size   = {[1 1] [1 -2]    [1 1]     [1 1]  [1 -1]     [1 1]    [1 1]   };

param = sw_readparam(inpForm, varargin{:});

nQ0     = param.nQ;
nWorker = param.nWorker;
nThread = param.nThread;
nRun    = param.nRun;
fName   = param.fName;
nSlice  = param.nSlice;
hermit  = param.hermit;

% setup
swpref.setpref('usemex',false,'tid',0,'fid',0);
Q = [40 40 nQ0];
yig = yig_create;

if nThread > 0
    setenv('OMP_NUM_THREADS',num2str(nThread));
else
    setenv('OMP_NUM_THREADS','');
end
% print header
measfun;

% runs without parallel pool
evalc('delete(gcp(''nocreate''))');
measfun(@eig_omp_duc, {yig Q 'hermit', hermit},false,nSlice,nRun,fName);
measfun(@eig_omp_duc, {yig Q 'hermit', hermit},true, nSlice,nRun,fName);

for ii = 1:numel(nWorker)
    nQ = round(nQ0/nWorker(ii))*nWorker(ii);
    Q  = rand(3,nQ);

    % run with parpool
    evalc(['parpool(' num2str(nWorker(ii)) ')']);
    measfun(@eig_omp_duc, {yig Q 'hermit', hermit},false, nSlice,nRun,fName);
    
    % stop pool
    evalc('delete(gcp(''nocreate''))');
end

end