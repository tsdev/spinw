function specOut = spinwave_spmd(obj,hkl,varargin)

% check parallel pool
pPool = gcp('nocreate');
if isempty(pPool)
    numWorker = 1;
else
    numWorker = pPool.NumWorkers;
end

Q      = sw_qscan(hkl);
nPoint = size(Q,2);
nPoint = round(nPoint/numWorker)*numWorker;
Q      = Q(:,1:nPoint);


Qc = Composite();
Ni = nPoint/numWorker;

for ii = 1:numWorker
    Qc{ii} = Q(:,(1:Ni)+(ii-1)*Ni);
end
spmd
    swpref.setpref('fid',0,'tid',0);
    spec = spinwave(obj,Qc,varargin{:});
    Sab  = spec.Sab;
    om   = spec.omega;
end

specOut       = spec{1};
specOut.Sab   = cat(4,Sab{:});
specOut.omega = cat(2,om{:});

end