function specOut = spinwavefast_duc_spmd(obj,hkl,varargin)

% check parallel pool
pPool = gcp('nocreate');
if isempty(pPool)
    numWorker = 1;
else
    numWorker = pPool.NumWorkers;
end

Q      = sw_qscan(hkl);
nPoint = size(Q,2);
nPoint = floor(nPoint/numWorker)*numWorker;
Q      = Q(:,1:nPoint);

% number of cores
nCore   = feature('numcores');
% number of threads per worker
nThread = floor(nCore*2/numWorker);

Qc = Composite();
Ni = nPoint/numWorker;

% data to transfer to workers
data = {obj varargin nThread};

for ii = 1:numWorker
    Qc{ii} = Q(:,(1:Ni)+(ii-1)*Ni);
end
spmd
    % set number of threads
    warning('off','MATLAB:maxNumCompThreads:Deprecated');
    maxNumCompThreads(data{3});
    swpref.setpref('fid',0,'tid',0);
    spec  = spinwavefast_duc(data{1},Qc,data{2}{:});
    Sperp = spec.Sperp;
    om    = spec.omega;
end

specOut       = spec{1};
specOut.Sperp = cat(4,Sperp{:});
specOut.omega = cat(2,om{:});

end