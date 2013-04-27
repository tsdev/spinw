function optm = optmagstr(obj, varargin)
% optimises magnetic structure by minimizing the energy using non-linear optimization algorithms.
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
% npar      Number of optimised parameters, optional. At least one of
%           {xmin,xmax,x0,npar} has to be defined.
% boundary  Boundary conditions of the extended unit cell.
%               'free'  Free, interactions between extedned unit cells are
%                       omitted.
%               'per'   Periodic, interactions between extended unit cells
%                       are retained.
%           Default is {'per' 'per' 'per'}.
% epsilon   The smalles value of incommensurability that is tolerated
%           without warning. Default is 1e-5.
%
% Optimisation parameters:
%
% tolx          Minimum change of x when convergence reached, default 
%               value is 1e-4.
% tolfun        Minimum change of the R value when convergence reached,
%               default value is 1e-5.
% maxfunevals   Maximum number of function evaluations, default value
%               is 1e7.
%
%
% Output:
%
% optm is struct type with the following fields:
% obj       sw object that contains the optimised magnetic structure.
% x         Optimised paramters.
% e         Energy per spin in the optimised structure.
% exitflag  Exit flag of the optimisation code, see fminsearch.
% output    Detailed output of the optimisation code, see fminsearch.
% param     Input parameters, stored in a struct.
%
% See also SW, SW.ANNEAL, GM_SPHERICAL3D, GM_PLANAR, FMINSEARCH.
%

inpForm.fname  = {'epsilon' 'func'           'boundary'          'xmin'   'xmax'  'x0'    'nPar'};
inpForm.defval = {1e-5      @gm_spherical3d  {'per' 'per' 'per'} []       []      []      []    };
inpForm.size   = {[1 1]     [1 1]            [1 3]               [1 -1]   [1 -2]  [1 -3]  [1 1] };
inpForm.soft   = {0         0                0                   1        1       1       1     };

inpForm.fname  = [inpForm.fname  {'tolx' 'tolfun' 'maxfunevals'}];
inpForm.defval = [inpForm.defval {1e-4   1e-5     1e7          }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]    [1 1]        }];
inpForm.soft   = [inpForm.soft   {0      0        0            }];

param = sw_readparam(inpForm, varargin{:});

[SS, SI] = obj.intmatrix;

if isempty(param.nPar)
    nPar = max(max(length(param.xmin),length(param.xmax)),length(param.x0));
else
    nPar = param.nPar;
end

% Initial parameters are random if param.x0 is undefined.
if ~isempty(param.x0)
    x0 = param.x0;
elseif isempty(param.xmax) && isempty(param.xmin)
    x0 = rand(1,nPar);
elseif isempty(param.xmax)
    x0 = param.xmin + rand(1,nPar);
elseif isempty(param.xmin)
    x0 = param.xmax - rand(1,nPar);
else
    x0 = rand(1,nPar).*(param.xmax-param.xmin)+param.xmin;
end
param.x0 = x0;

% creat initial magnetic structure
obj.genmagstr(param);




S       = sqrt(sum(obj.mag_str.S.^2,1));
nExt    = double(obj.mag_str.N_ext);
nMagExt = length(S);

% Modify the interaction matrices according to the boundary conditions.
for ii = 1:3
    if strcmp('free',param.boundary{ii})
        SS.all(:,SS.all(ii,:)~=0) = [];
    end
end

dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
atom1 = [SS.all(4,:)   1:nMagExt];
atom2 = [SS.all(5,:)   1:nMagExt];
JJ    = cat(3,reshape(SS.all(6:end,:),3,3,[]),SI.aniso);

Bgamma = obj.unit.gamma * SI.field;

[x, ~, exitflag, output] = sw_fminsearchbnd(@(x)efunc(x, S, dR, atom1, atom2, JJ, nExt, Bgamma, param.epsilon, param.func),x0,param.xmin,param.xmax,...
    optimset('TolX',param.tolx,'TolFun',param.tolfun,'MaxFunEvals',param.maxfunevals,'Display','off'));

[M, k, n] = param.func(S, x);

e = efunc(x, S, dR, atom1, atom2, JJ, nExt, Bgamma, param.epsilon, param.func);

obj.mag_str.S = M;
obj.mag_str.k = k;
obj.mag_str.n = n;

validate(obj);

% Create output struct
optm.obj      = copy(obj);
optm.x        = x;
optm.e        = e;
optm.exitflag = exitflag;
optm.output   = output;
optm.param    = param;

end

%% Energy function
function E = efunc(x, S, dR, atom1, atom2, JJ, nExt, Bgamma, epsilon, func)
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

E =  (sum(sum(sum(Ml.*JJ.*Mr))) - sum(Bgamma*M))/nMagExt;

end