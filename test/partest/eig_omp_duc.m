function spec = eig_omp_duc(~,M,varargin)

inpForm.fname  = {'hermit' 'optmem'};
inpForm.defval = {true     false};
inpForm.size   = {[1 1]   [1 1]};

param = sw_readparam(inpForm, varargin{:});
useMex = swpref.getpref('usemex',[]);

% check parallel pool
pPool = gcp('nocreate');
if isempty(pPool)
    numWorker = 0;
else
    numWorker = pPool.NumWorkers;
end

if numWorker==0
    
    M = num2cell(M);
    M = rand(M{:});
    if param.hermit
        M = (M+permute(M,[2 1 3]))/2;
    end
    spec = struct('V',{},'D',{});
    [spec.V,spec.D] = eigorth(M, 1e-5, true, useMex);
    
else
    M(3) = round(M(3)/numWorker)*numWorker;
    
    M = num2cell(M);
    M = rand(M{:});
    if param.hermit
        M = (M+permute(M,[2 1 3]))/2;
    end
    
    Mc = Composite();
    Mi = M(3)/numWorker;
    
    for ii = 1:numWorker
        Mc{ii} = M(:,:,(1:Mi)+(ii-1)*Mi);
    end
    spmd(nn)
        [V,D] = eigorth(Mc, tol, true, useMex);
    end
    
    spec   = struct;
    spec.V = cat(3,V{:});
    spec.D = cat(2,D{:});
end
end