function [pOpt,fVal,stat] = simplex(dat,func,p0,varargin)
% simplex optimisation
%
% ### Syntax
%
% `[pOpt,fVal,stat] = ndbase.simplex([],func,p0,Name,Value)`
%
% `[pOpt,fVal,stat] = ndbase.simplex(dat,func,p0,Name,Value)`
%
% ### Description
% 
% `[pOpt,fVal,stat] = ndbase.simplex([],func,p0,Name,Value)` finds a
% minimum of a function of several parameters using the constrained simplex
% optimization algorithm by calling the unconstrained Matlab built-in
% algorithm [matlab.fminsearch].
%
% The function has two different modes, depending on the first input
% argument. If `dat` is empty, `simplex` minimizes the cost function func,
% that has the following header:
% ```
% R2 = func(p)
% ```
%
% If `dat` is a struct, `simplex` optimizes the model defined by `func` via
% least squares to the data stored in `dat`. In this case `func` has the
% following header:
% ```
% yModel = func(x,p)
% ```
%
% And the least square deviation is defined by:
%
% $R^2 = \sum \frac{(y_{dat}-y_{fit})^2}{\sigma_{dat}^2}$
%
%  {{note If options is supplied, then TolX will apply to the transformed
%  variables. All other [matlab.fminsearch] parameters should be unaffected.
%
%  Variables which are constrained by both a lower and an upper
%  bound will use a $\sin(x)$ transformation. Those constrained by
%  only a lower or an upper bound will use a quadratic
%  transformation, and unconstrained variables will be left alone.
%
%  Variables may be fixed by setting their respective bounds equal.
%  In this case, the problem will be reduced in size for [matlab.fminsearch].
%
%  The bounds are inclusive inequalities, which admit the
%  boundary values themselves, but will not permit any function
%  evaluations outside the bounds. These constraints are strictly
%  followed.
%
%  If your problem has an exclusive (strict) constraint which will
%  not admit evaluation at the bound itself, then you must provide
%  a slightly offset bound. An example of this is a function which
%  contains the log of one of its parameters. If you constrain the
%  variable to have a lower bound of zero, then `simplex` may
%  try to evaluate the function exactly at zero.}}
%
% ### Examples
%
% Example usage on the rosen function.
%
% ```
% >>rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2
% ```
%
% Unconstrained optimisation:
%
% ```
% >>>rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2
% >>fminsearch(rosen,[3 3])>>
% ```
%
% Constrained optimisation:
%
% ```
% >>>rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2
% >>ndbase.simplex([],rosen,[3 3],'lb',[2 2],'ub',[])>>
% ```
%
% ### Input Arguments
%
% `dat`
% : Either empty or contains data to be fitted stored in a structure with
%   fields:
%   * `dat.x`   vector of $N$ independent variables,
%   * `dat.y`   vector of $N$ data values to be fitted,
%   * `dat.e`   vector of $N$ standard deviation (positive numbers)
%               used to weight the fit. If zero or missing
%               `1/dat.y^2` will be assigned to each point.
%
% `func`
% : Function handle with one of the following definition:
%   * `R2 = func(p)`        if `dat` is empty,
%   * `y  = func(x,p)`      if `dat` is a struct.
%   Here `x` is a vector of $N$ independent variables, `p` are the
%   $M$ parameters to be optimized and `y` is the simulated model, `R2`
%   is the value to minimize.
%
% ### Name-Value Pair Arguments
%
% `'lb'`
% : Vector with $N$ elements, lower boundary of the parameters. Default
%   value is -inf.
%
% `'ub'`
% : Vector with $N$ elements, upper boundary of the parameters. Default
%   value is inf.
%
% `'MaxIter'`
% : Maximum number of iterations, default value is $100M$.
%
% `'MaxFunEvals'`
% : Maximum number of function evaluations, default value is
%   $1000M$. NOT IMPLEMENTED!
%
% `'TolX'`
% : Convergence tolerance for parameters, default value is $10^{-3}$.
%
% `'Display'`
% : Level of information to print onto the Command Line during
%   optimisation. Possible values are `'off'`, `'notify'` and `'iter'`.
%   Default value is `'off'`.
%
% `'TolFun'`
% : Convergence tolerance on the `R2` value (return value of `func` or
%   the weighted least square deviation from data). Default value is
%   $10^{-3}$.
%
%
% ### See Also
%
% [ndbase.lm] \| [ndbase.pso]

% Original author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 4
% Release date: 7/23/06

if nargin == 0
    swhelp ndbase.simplex
    return
end

% number of parameters
Np    = numel(p0);

% not implemented yet: 'MaxFunEvals'
inpForm.fname  = {'Display' 'TolFun' 'TolX' 'MaxIter' 'lb'      'ub'    };
inpForm.defval = {'off'     1e-3     1e-3   100*Np    []        []      };
inpForm.size   = {[1 -1]    [1 1]    [1 1]  [1 1]     [-5 -2]   [-3 -4] };
inpForm.soft   = {false     false    false  false     true      true    };

param = sw_readparam(inpForm, varargin{:});
param.Np = Np;

% limits
LB = param.lb;
UB = param.ub;

if ~isempty(LB) && ~isempty(UB) && any(UB<LB)
    error('simplex:WrongInput','Upper boundary has to be larger than the lower boundary!');
end

% number of free parameters
if isempty(LB) || isempty(UB)
    Nv = max(numel(LB),numel(UB));
else
    Nv = sum(UB>LB);
end

% check input function
if ischar(func)
    % convert to fuction handle
    func = str2func(func);
end

if ~isa(func,'function_handle')
    error('simplex:WrongInput','The input function is neither a string, not function handle!');
end

% define weighted least squares if dat is given
if ~isempty(dat)
    dat.x = dat.x(:);
    dat.y = dat.y(:);
    
    if ~isfield(dat,'e') || isempty(dat.e) || ~any(dat.e(:))
        weight = 1./abs(dat.y);
    else
        if any(dat.e(:)<0)
            error('pso:WrongInput','Standard deviations have to be positive!')
        end
        weight = 1./dat.e(:).^2;
    end
    func0 = func;
    func = @(p)sum(weight.*(func(dat.x(:),p)-dat.y(:)).^2);
end

if isempty(LB)
  LB = repmat(-inf,Np,1);
else
  LB = LB(:);
end
if (nargin<4) || isempty(UB)
  UB = inf(Np,1);
else
  UB = UB(:);
end

if (Np~=numel(LB)) || (Np~=numel(UB))
  error('simplex:WrongInput','p0 is incompatible in size with the given limits!')
end


% 0 --> unconstrained variable
% 1 --> lower bound only
% 2 --> upper bound only
% 3 --> dual finite bounds
% 4 --> fixed variable
param.BoundClass = zeros(Np,1);
for ii=1:Np
  k = isfinite(LB(ii)) + 2*isfinite(UB(ii));
  param.BoundClass(ii) = k;
  if (k==3) && (LB(ii)==UB(ii))
    param.BoundClass(ii) = 4;
  end
end

% transform starting values into their unconstrained
% surrogates. Check for infeasible starting guesses.
p0u = p0;
k   = 1;
for ii = 1:Np
  switch param.BoundClass(ii)
    case 1
      % lower bound only
      if p0(ii)<=LB(ii)
        % infeasible starting value. Use bound.
        p0u(k) = 0;
      else
        p0u(k) = sqrt(p0(ii) - LB(ii));
      end
      
      % increment k
      k=k+1;
    case 2
      % upper bound only
      if p0(ii)>=UB(ii)
        % infeasible starting value. use bound.
        p0u(k) = 0;
      else
        p0u(k) = sqrt(UB(ii) - p0(ii));
      end
      
      % increment k
      k=k+1;
    case 3
      % lower and upper bounds
      if p0(ii)<=LB(ii)
        % infeasible starting value
        p0u(k) = -pi/2;
      elseif p0(ii)>=UB(ii)
        % infeasible starting value
        p0u(k) = pi/2;
      else
        p0u(k) = 2*(p0(ii) - LB(ii))/(UB(ii)-LB(ii)) - 1;
        % shift by 2*pi to avoid problems at zero in fminsearch
        % otherwise, the initial simplex is vanishingly small
        p0u(k) = 2*pi+asin(max(-1,min(1,p0u(k))));
      end
      
      % increment k
      k=k+1;
    case 0
      % unconstrained variable. x0u(i) is set.
      p0u(k) = p0(ii);
      
      % increment k
      k=k+1;
    case 4
      % fixed variable. drop it before fminsearch sees it.
      % k is not incremented for this variable.
  end
  
end
% if any of the unknowns were fixed, then we need to shorten
% x0u now.
if k<=Np
  p0u(k:Np) = [];
end

% were all the variables fixed?
if isempty(p0u)
  % All variables were fixed. quit immediately, setting the
  % appropriate parameters, then return.
  
  % undo the variable transformations into the original space
  x = xtransform(p0u,param);
  
  % stuff fval with the final value
  fVal = func(x,p0);
  
  stat            = struct;
  stat.msg        = 'Parameters are fixed, no optimisation';
  stat.p          = p0;
  stat.sigP       = [];
  stat.Rsq        = [];
  stat.sigY       = [];
  stat.corrP      = [];
  stat.cvgHst     = [];
  stat.iterations = 0;
  stat.funcCount  = 1;
  stat.algorithm  = 'Nelder-Mead simplex direct search';
  stat.exitFlag   = 0;
  stat.param = param;
  
  if isempty(dat)
      stat.func   = func;
  else
      stat.func   = func0;
  end
  % return with no call at all to fminsearch
  return
end


% now we can call fminsearch, but with our own
% intra-objective function.
intrafun = @(x)func(xtransform(x,param));

[pu,fVal,exitFlag,stat0] = fminsearch(intrafun,p0u,param);

% undo the variable transformations into the original space
pOpt = xtransform(pu,param);

stat            = struct;
stat.p          = pOpt;
stat.sigP       = [];
if isempty(dat)
    stat.redX2 = fVal;
else
    % divide R2 with the statistical degrees of freedom
    stat.redX2   = fVal/(numel(dat.x)-Nv+1);
end
stat.msg        = stat0.message;
stat.Rsq        = [];
stat.sigY       = [];
stat.corrP      = [];
stat.cvgHst     = [];
stat.iterations = stat0.iterations;
stat.funcCount  = stat0.funcCount;
stat.algorithm  = stat0.algorithm;
stat.exitFlag   = exitFlag;
stat.param = param;

if isempty(dat)
    stat.func   = func;
else
    stat.func   = func0;
end

end % mainline end

% ======================================
function xtrans = xtransform(x,p)
% converts unconstrained variables into their original domains

xtrans = x*0;
% k allows some variables to be fixed, thus dropped from the
% optimization.
k=1;
for i = 1:p.Np
  switch p.BoundClass(i)
    case 1
      % lower bound only
      xtrans(i) = p.lb(i) + x(k).^2;
      
      k=k+1;
    case 2
      % upper bound only
      xtrans(i) = p.ub(i) - x(k).^2;
      
      k=k+1;
    case 3
      % lower and upper bounds
      xtrans(i) = (sin(x(k))+1)/2;
      xtrans(i) = xtrans(i)*(p.ub(i) - p.lb(i)) + p.lb(i);
      % just in case of any floating point problems
      xtrans(i) = max(p.lb(i),min(p.ub(i),xtrans(i)));
      
      k=k+1;
    case 4
      % fixed variable, bounds are equal, set it at either bound
      xtrans(i) = p.lb(i);
    case 0
      % unconstrained variable.
      xtrans(i) = x(k);
      
      k=k+1;
  end
end

end % sub function xtransform end