function varargout = optmagsteep2(obj, varargin)
% optimise magnetic structure using the steepest descendend method
%
% OPTMAGSTEEP(obj, 'option1', value1 ...)
%
% The function cannot deal with single ion anisotropy and incommensurate
% structures at the present state. These features will come soon.
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
%
% Output:
%
% E         If requested, calculates the energy after every iteration, the
%           calculation makes the script very slow.
%
%
% See also SW, SW.ANNEAL, SW_FSUB, SW_FSTAT.
%

nExt   = double(obj.mag_str.N_ext);

inpForm.fname  = {'nRun' 'epsilon' 'random' 'boundary'          'subLat'};
inpForm.defval = {100     1e-5      false   {'per' 'per' 'per'}  []     };
inpForm.size   = {[1 1]   [1 1]     [1 1]   [1 3]                [1 -1] };
inpForm.soft   = {0       0         0       0                    1      };

inpForm.fname  = [inpForm.fname  {'nExt' 'fSub'   }];
inpForm.defval = [inpForm.defval {nExt   @sw_fsub }];
inpForm.size   = [inpForm.size   {[1 3]  [1 1]    }];
inpForm.soft   = [inpForm.soft   {0      0        }];


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
M0  = obj.mag_str.S;

% Produce the interaction matrices
[SS, SI] = obj.intmatrix;

% Function options.
nRun    = param.nRun;
nMagExt = size(M0,2);

% Initial moment directions, the extra zeros are for easy indexing, size of
% M is (1,nMagExt+1).
M     = [M0 [0;0;0]];
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

% Number of moments on each sublattice.
nElementSub = sum(Sindex,2);

% Speeds up the code by storing every sublattice data in different cells
%csSSJ  = cell(nSub,1);
%csSSi  = cell(nSub,1);
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
    E = zeros(1,nRun);
end

for ii = 1:nRun
    for jsub = 1:nSub
        % Logical vector, selecting the moments on a given
        % sublattice [1,nMagExt]
        sSindex = Sindex(jsub,:);
        % Stores the interaction strength
        % SSJ(:,sSindex) [maxNeigh, nElementSub(jsub)]
        % --> reshape into [1, maxNeigh*nElementSub(jsub)] size
        % --> repmat into nSpinDim vertically
        % --> final size [nSpinDim, maxNeigh*nElementSub(jsub)]
        %sSSJ = csSSJ{jsub};
        % Stores the indices of neighbouring magnetic moments.
        % SSi(:,sSindex) [maxNeigh, nElementSub(jsub)]
        % --> reshape into [1, maxNeigh*nElementSub(jsub)] size
        %sSSi = csSSi{jsub};
        % Previous moment directions of the given sublattice.
        %MOld = M(:,sSindex);
        % F stores the molecular field acting on the moments of
        % the jsub sublattice (exchange+external field).
        % F [3, nElementSub(jsub)]
        sSSJG = csSSJG{jsub};
        F = squeeze(sum(reshape(permute(mmat(sSSJG,permute(M(:,csSSiG{jsub}),[1 3 2])),[1 3 2]),3, maxNeighG,[]),2));
        
        % Adds external magnetic field.
        if param.isfield
            F = F - cB{jsub};
        end
        
        % Generate new spin directions, creating normal
        % distribution of coordinates, then normalizing them.
        if ~param.aniso
            MNew  = -bsxfun(@times,F,cS{jsub}./sqrt(sum(F.^2)));
        else
            % Calculates old anisotropy field.
            %AOld = [sum(MOld.*cAx{jsub},1); sum(MOld.*cAy{jsub},1); sum(MOld.*cAz{jsub},1)];
            % Calculates new anisotropy field.
            %ANew = [sum(MNew.*cAx{jsub},1); sum(MNew.*cAy{jsub},1); sum(MNew.*cAz{jsub},1)];
            % Calculates the energy difference on each spin.
            %dE = sum((MNew-MOld).*F+MNew.*ANew-MOld.*AOld,1);
        end
        
        aidx = true(1,nElementSub(jsub));
        % sidx stores the accepted spin indices in the whole spin list.
        sidx = false(nMagExt+1,1);
        sidx(sSindex) = aidx;
        
        M(:,sidx) = MNew(:,aidx);
        
    end
    
    
    % Calculates the system energy at the end of the temperature step.
    obj.mag_str.S = M(:,1:end-1);
    if nargout == 1
        E(ii) = obj.energy;
    end
    
    if fid == 1
        sw_status(ii/param.nRun*100);
    end
    
end

if fid == 1
    sw_status(100,2);
else
    if fid ~= 0
        fprintf0(fid,'Calculation finished.\n');
    end
end

obj.mag_str.S = M(:,1:end-1);

if nargout == 1
    varargout{1} = E;
end

end