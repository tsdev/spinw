function spec = eig_omp_duc(M,varargin)

inpForm.fname  = {'hermit'};
inpForm.defval = {true    };
inpForm.size   = {[1 1]   };

param = sw_readparam(inpForm, varargin{:});
useMex = swpref.getpref('usemex',[]);

% check parallel pool
pPool = gcp('nocreate');
if isempty(pPool)
    numWorker = 1;
else
    numWorker = pPool.NumWorkers;
end

M(3) = round(M(3)/numWorker)*numWorker;

M = num2cell(M);
M = rand(M{:});
if param.hermit
    M = (M+permute(M,[1 3 2]))/2;
end

% [spec.V,spec.D] = eigorth(M, tol, true, useMex);
% 
% spec.param.nSlice = 1;
% spec.param.hermit = param.hermit;


Mc = Composite();
Mi = M(3)/numWorker;

for ii = 1:numWorker
    Mc{ii} = M(:,:,(1:Mi)+(ii-1)*Mi);
end
spmd
    [V,D] = eigorth(Mc, tol, true, useMex);
end

spec   = struct;
spec.V = cat(3,V{:});
spec.D = cat(2,D{:});

end