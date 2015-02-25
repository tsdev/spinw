function optm = optmagsteep(obj, varargin)
% optimise magnetic structure using the steepest descendend method
%
% optm = OPTMAGSTEEP(obj, 'option1', value1 ...)
%
% The function works only for magnetic structures with zero k-vector.
%
%
% Input:
%
% obj             Input object contains structural data, sw type.
%
% Options:
%
% nRun      Number of iterations, default is 100 (it is usually enough).
% boundary  Boundary conditions of the extended unit cell.
%               'free'  Free, interactions between extedned unit cells are
%                       omitted.
%               'per'   Periodic, interactions between extended unit cells
%                       are retained.
%           Default is {'per' 'per' 'per'}.
% nExt      The size of the magnetic cell in number of unit cells, to
%           provide input information to 'fStat'.
%           Default is from obj.mag_str.N_ext.
% fSub      Function to define sublattices for Monte-Carlo speedup.
%           cGraph = fSub(conn,nExt), where cGraph is a (1,nMagExt) sized
%           vector, conn is a (2,nConn) size matrix and nExt is equal to
%           'nExt'. Default is <a href="matlab: doc sw_fsub">@sw_fsub</a>
% subLat    Vector that assigns all magnetic moments into non-interacting
%           sublattices, contains a single index (1,2,3...) for every
%           magnetic moment, size is (1,nMagExt). If undefined, the
%           function defined in 'fSub' will be used to partition the
%           lattice.
% random    Random initial conditions, if initial spin configuration
%           is undefined (obj.mag_str.S is empty) the initial configuration
%           is automaticly random independently of the value of random.
%           Default is false.
% TolX      Minimum change of the magnetic moment when the algorithm stops.
% saveAll   Save moment directions for every loop, default is false.
%
% Output:
%
% 'optm' is a struct type variable with the following fields:
% obj       sw object that contains the optimised magnetic structure.
% M         Magnetic moment directions with dimensions [3 nMagExt], if
%           'saveAll' parameter is true, it contains the magnetic structure
%           after every loop in a matrix with dimensions [3 nMagExt nLoop].
% dM     	The change of magnetic moment vector averaged over all moments
%           in the last loop.
% e         Energy per spin in the optimised structure.
% param     Input parameters, stored in a struct.
% nRun      Number of loops executed.
% datestart Starting time of the function.
% dateend   End time of the function.
% title     Title of the simulation, given in the input.
%
% See also SW, SW.ANNEAL, SW_FSUB, SW_FSTAT.
%

% save the time of the beginning of the calculation
if nargout > 0
    optm.datestart  = datestr(now);
end

nExt   = double(obj.mag_str.N_ext);
title0 = 'Optimised magnetic structure using steepest descendent';

inpForm.fname  = {'nRun' 'epsilon' 'random' 'boundary'          'subLat'};
inpForm.defval = {100     1e-5      false   {'per' 'per' 'per'}  []     };
inpForm.size   = {[1 1]   [1 1]     [1 1]   [1 3]                [1 -1] };
inpForm.soft   = {0       0         0       0                    1      };

inpForm.fname  = [inpForm.fname  {'nExt' 'fSub'   'TolX' 'title' 'saveAll'}];
inpForm.defval = [inpForm.defval {nExt   @sw_fsub 1e-10  title0  false    }];
inpForm.size   = [inpForm.size   {[1 3]  [1 1]    [1 1]  [1 -2]  [1 1]    }];
inpForm.soft   = [inpForm.soft   {0      0        0      0        0        }];


param = sw_readparam(inpForm,varargin{:});

% Text output file
fid = obj.fid;

% Creates random spin directions if param.random is true.
mag_param = struct;
if param.random
    mag_param.mode = 'random';
else
    mag_param.mode = 'extend';
end
mag_param.nExt = param.nExt;

obj.genmagstr(mag_param);
M  = obj.mag_str.S;

% Produce the interaction matrices
[SS, SI] = obj.intmatrix;

% Function options.
nRun    = param.nRun;
nMagExt = size(M,2);

% Spin length for normalization.
S     = sqrt(sum(M.^2,1));

% Modify the interaction matrices according to the boundary conditions.
for ii = 1:3
    if strcmp('free',param.boundary{ii})
        SS.all(:,SS.all(ii,:)~=0) = [];
    end
end

% Since k_m=(0,0,0) the spins that are coupled to themself contribute with
% a constant self-energy, removing this doesn't change thermodynamical
% behaviour just shifts the zero energy.
SS.all(:, SS.all(4,:)==SS.all(5,:)) = [];

% Calculates the energy of the initial configuration and prepares the
% anisotropy matrix. B is in units of the couplings.
Bloc = permute(mmat(SI.field*obj.unit.muB,SI.g),[2 3 1]);
AA = SI.aniso;
Ax = squeeze(AA(:,1,:));
Ay = squeeze(AA(:,2,:));
Az = squeeze(AA(:,3,:));

% Checks whether there is any external field
param.isfield = any(Bloc(:));

% Checks whether anisotropy is non-zero.
if any(AA(:))
    param.aniso = true;
else
    param.aniso = false;
end

% Assing moments to sublattices for parallel calculation. There are no
% coupling between moments on the same sublattice, thus Weiss field can be
% calculated parallel. SSc stores the index of the sublattice, size:
% (1,nMagExt)

if isempty(param.subLat)
    SSc  = param.fSub(SS.all(4:5,:),param.nExt);
    param.subLat = SSc;
else
    SSc  = param.subLat;
end
nSub = max(SSc);

nNeighG = zeros(nMagExt,1);
for ii = 1:nMagExt
    nNeighG(ii) = sum((SS.all(4,:) == ii)|(SS.all(5,:) == ii));
end

% Maximum number of neighbours
maxNeighG = max(nNeighG);

% Interaction matrices and neigbor indices
SSiG = zeros(maxNeighG,nMagExt) + (nMagExt+1);
SSJG = zeros(9,maxNeighG,nMagExt);

% Indexes for transposing J for exchanged spins in the interaction.
% Default is Si * J * Sj, or Sj * J' * Si has to be used.
trIdx = reshape(reshape(1:9,[3 3])',[9 1])+5;
for ii = 1:nMagExt
    SSiG(1:nNeighG(ii),ii)   = [SS.all(5,(SS.all(4,:) == ii))    SS.all(4,(SS.all(5,:) == ii))    ]';
    SSJG(:,1:nNeighG(ii),ii) = [SS.all(6:14,(SS.all(4,:) == ii)) SS.all(trIdx,(SS.all(5,:) == ii))];
end

% Store spin indices of each sublattice for speedup.
Sindex = zeros(nSub,nMagExt);
Sindex(nSub*(0:nMagExt-1)+SSc) = 1;
Sindex      = logical(Sindex);

% Speeds up the code by storing every sublattice data in different cells
csSSiG = cell(nSub,1);
csSSJG = cell(nSub,1);
cAx    = cell(nSub,1);
cAy    = cell(nSub,1);
cAz    = cell(nSub,1);
cS     = cell(nSub,1);
cB     = cell(nSub,1);

for ii = 1:nSub
    sSindex    = Sindex(ii,:);
    cS{ii}    = S(sSindex);
    
    cAx{ii}   = Ax(:,sSindex);
    cAy{ii}   = Ay(:,sSindex);
    cAz{ii}   = Az(:,sSindex);
    
    csSSiG{ii} = reshape(SSiG(:,sSindex),1,[]);
    csSSJG{ii} = reshape(SSJG(:,:,sSindex),3,3,[]);
    
    cB{ii}    = Bloc(:,sSindex);
end

if fid == 1
    sw_status(0,1);
end

if nargout == 1
    E   = zeros(1,nRun);
    if param.saveAll
        Msave = zeros(3,nMagExt);
    end
end

% Initial step size is infinite, to make at least 1 cycle.
dM = inf;
% Initial index.
rIdx = 0;

while (rIdx < nRun) && (dM>param.TolX)
    Mold = M;
    for jsub = 1:nSub
        % Logical vector, selecting the moments on a given
        % sublattice [1,nMagExt]
        sSindex = Sindex(jsub,:);
        % F stores the molecular field acting on the moments of
        % the jsub sublattice (exchange+external field).
        % F [3, nElementSub(jsub)]
        sSSJG = csSSJG{jsub};
        F = squeeze(sum(reshape(permute(mmat(sSSJG,permute(M(:,csSSiG{jsub}),[1 3 2])),[1 3 2]),3, maxNeighG,[]),2));
        
        % Adds external magnetic field.
        if param.isfield
            F = F - cB{jsub};
        end
        
        % Adds anisotropy field.
        if param.aniso
            % Select the moment vectors on the sublattice.
            Ms = M(:,sSindex);
            Fa = 2*[sum(Ms.*cAx{jsub},1); sum(Ms.*cAy{jsub},1); sum(Ms.*cAz{jsub},1)];
            F = F + Fa;
        end
        
        M(:,sSindex) = -bsxfun(@times,F,cS{jsub}./sqrt(sum(F.^2)));
        
    end
    
    % Calculates the system energy at the end of the temperature step.
    if nargout > 0
        obj.mag_str.S = M;
        E(rIdx+1) = obj.energy;
        if param.saveAll
            Msave(:,:,rIdx+1) = M;
        end
    end
    
    if fid == 1
        sw_status(rIdx/param.nRun*100);
    end
    
    % Check stopping condition, give the dM limit.
    dM = sum(sqrt(sum((Mold - M).^2,1)))/nMagExt;
    rIdx = rIdx + 1;
    
end

if fid == 1
    sw_status(100,2);
else
    if fid ~= 0
        fprintf0(fid,'Calculation finished.\n');
    end
end

% Save optimised magnetic structure into the sw object.
obj.mag_str.S = M;

% Create output structure.
if nargout > 0
    optm.obj      = copy(obj);
    if param.saveAll
        optm.M = Msave;
    else
        optm.M = M;
    end
    optm.dM       = dM;
    optm.e        = E(1:rIdx);
    optm.param    = param;
    optm.nRun     = rIdx;
    optm.dateend  = datestr(now);
    optm.title    = param.title;
end

end