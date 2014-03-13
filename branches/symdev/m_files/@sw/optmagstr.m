function optm = optmagstr(obj, varargin)
% optimises magnetic structure by minimizing the energy using non-linear optimization algorithms
%
% optm = OPTMAGSTR(swobj, Option1, Value1, ...)
%
% Options:
%
% func      Function that produce the magnetic moments, ordering wave
%           vector and normal vector from the optimization
%           parameters in the following form:
%               [M, k, n] = @(x)func(M0, x)
%           where M is (3,nMagExt) size matrix. k is the ordering
%           wave vector, its size is (1,3). n is the normal vector
%           of the spin rotation plane, its size is (1,3). The
%           default is @gm_spherical3d. For planar magnetic structure
%           use @gm_planar.
% xmin      Minimum limit of the optimisation parameters, optional.
% xmax      Maximum limit of the optimisation parameters, optional.
% x0        Starting value of the optimisation parameters. If empty
%           or undefined, then random values are used.
% boundary  Boundary conditions of the extended unit cell.
%               'free'  Free, interactions between extedned unit cells are
%                       omitted.
%               'per'   Periodic, interactions between extended unit cells
%                       are retained.
%           Default is {'per' 'per' 'per'}.
% epsilon   The smalles value of incommensurability that is tolerated
%           without warning. Default is 1e-5.
% nRun      Number of runs. If random starting parameters are given, the
%           optimisation process will be rerun nRun times and the best
%           result (lowest ground state energy per spin) will be saved in
%           the result.
% fid       For text output of the calculation:
%               0   No text output.
%               1   Output written onto the Command Window. (default)
%               fid Output written into a text file opened with the
%                   fid = fopen(path) command.
%
% Limits only on selected prameters:
%
% Limits can be given on any input parameter of the constraint function by
% giving the name of the parameter, see the help of the used constraint
% function in the following format: optmagstr('ParName',[min max],...).
% For example to fix the nTheta value of @gm_planar during the optimisation
% to zero use:
% optmagstr(obj,'func',@gm_planar,'nTheta',[0 0]);
%
% Optimisation parameters:
%
% tolx          Minimum change of x when convergence reached, default
%               value is 1e-4.
% tolfun        Minimum change of the R value when convergence reached,
%               default value is 1e-5.
% maxfunevals   Maximum number of function evaluations, default value
%               is 1e7.
% maxiter       Maximum number of iterations, default value is 1e4.
%
%
% Output:
%
% optm is struct type with the following fields:
% obj       sw object that contains the optimised magnetic structure.
% x         Optimised paramters, dimensions are [1 nPar].
% fname     Name of the contraint function.
% xname     Cell containing the name of the x parameters, dimensions are
%           [1 nPar].
% e         Energy per spin in the optimised structure.
% exitflag  Exit flag of the optimisation code, see fminsearch.
% output    Detailed output of the optimisation code, see fminsearch.
% param     Input parameters, stored in a struct.
%
% See also SW, SW.ANNEAL, GM_SPHERICAL3D, GM_PLANAR, FMINSEARCH.
%

inpForm.fname  = {'epsilon' 'func'           'boundary'          'xmin'   'xmax'  'x0'   };
inpForm.defval = {1e-5      @gm_spherical3d  {'per' 'per' 'per'} []       []      []     };
inpForm.size   = {[1 1]     [1 1]            [1 3]               [1 -1]   [1 -2]  [1 -3] };
inpForm.soft   = {0         0                0                   1        1       1      };

inpForm.fname  = [inpForm.fname  {'tolx' 'tolfun' 'maxfunevals' 'nRun' 'maxiter' 'fid'}];
inpForm.defval = [inpForm.defval {1e-4   1e-5     1e7           1      1e4       1    }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]    [1 1]         [1 1]  [1 1]     [1 1]}];
inpForm.soft   = [inpForm.soft   {0      0        0             0      0         0    }];

param = sw_readparam(inpForm, varargin{:});

% creat initial magnetic structure
param.showWarn = false;
obj.genmagstr(param);

% starting magnetic structure from sw object
S       = sqrt(sum(obj.mag_str.S.^2,1));
nExt    = double(obj.mag_str.N_ext);
nMagExt = length(S);

% determine the limits from the constraint function
if nargout(param.func) == 6
    [~, ~, ~, fname, pname, limit] = param.func(S,[]);
    % limits of the fitting parameters
    nPar = size(limit,2);
    % limits of the parameters
    if numel(param.xmin) ~= nPar
        param.xmin = limit(1,:);
    end
    
    if numel(param.xmax) ~= nPar
        param.xmax = limit(2,:);
    end
    
    % check if any parameter name is given as explicit input with limits
    inpForm.fname  = pname;
    inpForm.defval = repmat({[]},1,nPar);
    inpForm.size   = repmat({[1 2]},1,nPar);
    inpForm.soft   = repmat({1},1,nPar);
    % test input parameters
    if nargin == 2
        varargin{1}.showWarn = false;
    else
        varargin{end+1} = 'showWarn';
        varargin{end+1} = false;
    end
    fparam = sw_readparam(inpForm, varargin{:});
    
    for ii = 1:nPar
        if ~isempty(fparam.(inpForm.fname{ii}))
            param.xmin(ii) = fparam.(inpForm.fname{ii})(1);
            param.xmax(ii) = fparam.(inpForm.fname{ii})(2);
        end
    end
    % restrict input parameter limits to the hardwired constraint function
    % limits
    param.xmin = max(param.xmin,limit(1,:));
    param.xmax = min(param.xmax,limit(2,:));
    
else
    if isempty(param.xmin) || isempty(param.xmax) || (numel(param.xmin) ~= numel(param.xmax))
        error('sw:optmagtr:WrongInput','Missing limits on the x fitting parameters (use xmin and xmax options)!');
    end
    nPar  = numel(param.xmin);
    fname = '';
    pname = repmat({''},1,nPar);
end

% get magnetic couplings
[SS, SI] = obj.intmatrix;

% Initial parameters are random if param.x0 is undefined/wrong size
xRand = (numel(param.x0)~= nPar);

if (~xRand) || (param.nRun<1)
    param.nRun = 1;
end

% Modify the interaction matrices according to the boundary conditions.
for ii = 1:3
    if strcmp('free',param.boundary{ii})
        SS.all(:,SS.all(ii,:)~=0) = [];
    end
end

dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
atom1 = [SS.all(4,:)   1:nMagExt];
atom2 = [SS.all(5,:)   1:nMagExt];
JJ    = cat(3,reshape(SS.all(6:14,:),3,3,[]),SI.aniso);
g     = SI.g;

% B * g: field * g-tensor
Bg  = permute(mmat(SI.field,g)*obj.unit.muB,[2 3 1]);

minE = 0;
minX = zeros(1,nPar);

fid = param.fid;

if fid == 1
    sw_status(0,1);
end

% Loop over nRun times
for ii = 1:param.nRun
    if xRand
        x0 = rand(1,nPar).*(param.xmax-param.xmin)+param.xmin;
    else
        x0 = param.x0;
    end
    [X, E, exitflag, output] = sw_fminsearchbnd(@(x)efunc(x, S, dR, atom1, atom2, JJ, nExt, Bg, param.epsilon, param.func),x0,param.xmin,param.xmax,...
        optimset('TolX',param.tolx,'TolFun',param.tolfun,'MaxFunEvals',param.maxfunevals,'MaxIter',param.maxiter,'Display','off'));
    if E < minE
        minE = E;
        minX = X;
    end
    if fid == 1
        sw_status(ii/param.nRun*100);
    end

end

if fid == 1
    sw_status(100,2);
else
    if fid ~= 0
        fprintf(fid,'Calculation finished.\n');
    end
end

[M, k, n] = param.func(S, minX);

obj.mag_str.S = M;
obj.mag_str.k = k;
obj.mag_str.n = n;

validate(obj);

% Create output struct
if nargout > 0
    optm.obj      = copy(obj);
    optm.x        = minX;
    optm.e        = minE;
    optm.exitflag = exitflag;
    optm.output   = output;
    optm.param    = param;
    optm.fname    = fname;
    optm.xname    = pname;
end

end

%% Energy function
function E = efunc(x, S, dR, atom1, atom2, JJ, nExt, Bg, epsilon, func)
% Cost function to optimize.

[M, k, n] = func(S, x);

M1 = M(:,atom1);
M2 = M(:,atom2);

kExt    = k.*nExt;
nMagExt = length(S);

% Incommenasurate structure in the extended unit cell.
incomm = any(abs(k-round(k))>epsilon);

% Rotate spins for incommensurate structures.
if incomm && any(any(dR))
    dRIdx = find(any(dR));
    for ii = 1:length(dRIdx)
        M2(:,dRIdx(ii)) = sw_rot(n, kExt*dR(:,dRIdx(ii))*2*pi, M2(:,dRIdx(ii)));
    end
end

Ml = repmat(permute(M1,[1 3 2]),[1 3 1]);
Mr = repmat(permute(M2,[3 1 2]),[3 1 1]);

E =  (sum(sum(sum(Ml.*JJ.*Mr))) - sum(Bg(:).*M(:)))/nMagExt;

end