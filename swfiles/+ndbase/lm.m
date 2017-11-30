function [p0,yHat,stat] = lm(dat,func0,p0,varargin)
% least square refinement of parameters using Levenberg Marquardt method
%
% [pOpt,fVal,stat] = NDBASE.LM(dat,func,p0,'Option1','Value1',...)
%
% Levenberg Marquardt curve-fitting function, minimizes sum of weighted
% squared residuals.
%
% Input:
%
% dat       Data to be fitted stored in a structure with fields:
%               dat.x   vector of N independent variables,
%               dat.y   vector of N data values to be fitted,
%               dat.e   vector of N standard deviation (positive numbers)
%                       used to weight the fit. If zero or missing
%                       1./dat.y^2 will be assigned to each point, however
%                       in this case the refinement is reduced to an
%                       optimisation problem.
% func      Function handle with the following definition:
%               y = func(x,p)
%           where x is a vector of N independent variables, p are the
%           M parameters to be optimized and y is the simulated model.
% p0        Vector of M initial parameters.
%
% Options:
%
% Options can be given using the modified output of optimset() or as option
% name string option value pairs.
%
% dp        Vector with N or 1 element, defines the fractional increment of
%           'p' when calculating the Jacobian matrix dFunc(x,p)/dp:
%               dp(j)>0     central differences calculated,
%               dp(j)<0     one sided 'backwards' differences calculated,
%               dp(j)=0     sets corresponding partials to zero, i.e. holds
%                           p(j) fixed.
%           Default value if 1e-3.
% vary      Vector with N elements, if an element is false, the
%           corresponding parameter will be fixed. Default value is 
%           false(1,N).
% win       Limits for the independent variabel values where the function
%           is fitted. Default is [-inf inf].
% lb        Vector with N elements, lower boundary of the parameters.
%           Default value is -inf.
% ub        Vector with N elements, upper boundary of the parameters.
%           Default value is inf.
% MaxIter   Maximum number of iterations, default value is 10*M.
% MaxFunEvals Maximum number of function evaluations, default value is
%           100*M.
% TolX      Convergence tolerance for parameters, defines the maximum of 
%           the relative chande of any parameter value. Default value is
%           1e-3.
% eps1      Convergence tolerance for gradient, default value is 1e-3.
% eps2      Convergence tolerance for reduced Chi-square, default value is
%           1e-2.
% eps3      Determines acceptance of a L-M step, default value is 0.1.
% lambda0   Initial value of L-M paramter, default value is 1e-2.
% nu0       Value that determines the speed of convergence. Default value
%           is 10. It should be larger than 1.
% lUp       Factor for increasing lambda, default value is 30.
% lDown     Factor for decreasing lambda, default value is 7.
% update    Type of parameter update:
%                   'lm'        Levenberg-Marquardt lambda update,
%                   'quadratic' Quadratic update,
%                   'nielsen'   Nielsen's lambda update equations (default).
% extraStat Calculates extra statistics: covariance matrix of parameters,
%           cofficient of multiple determination, asymptotic standard
%           error of the curve-fit and convergence history.
% confLev   Confidence level, where the error of the curve fit (sigY) is
%           calculated. Default is erf(1/sqrt(2))~0.6827 for standard
%           deviation (+/- 1 sigma).
%
% Output:
%
% pOpt      Value of the M optimal parameters.
% fVal      Value of the model function calculated with the optimal
%           parameters at the N independent values of x.
%
% stat      Structure, storing the detailed output of the calculation with
%           the following fields:
%               p       Least-squares optimal estimate of the parameter
%                       values.
%               redX2   Reduced Chi squared error criteria, its value
%                       should be close to 1. If the value is larger, the
%                       model is not a good description of the data. If the
%                       value is smaller, the model is overparameterized
%                       and fitting the statistical error of the data.
%               sigP    Asymptotic standard error of the parameters.
%               sigY    Asymptotic standard error of the curve-fit.
%               corrP   Correlation matrix of the parameters.
%               Rsq     R-squared cofficient of multiple determination.
%               cvgHst  Convergence history.
%               exitFlag The reason, why the code stopped:
%                           1       convergence in r.h.s. ("JtWdy"),
%                           2       convergence in parameters,
%                           3       convergence in reduced Chi-square,
%                           4       maximum Number of iterations reached
%                                   without convergence
%               message String, one of the above messages.
%               nIter   The number of iterations executed during the fit.
%               nFunEvals The number of function evaluations executed
%                       during the fit.
%
% See also NDBASE.PSO.
%

%   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. 4 May 2016
%   modified from: http://octave.sourceforge.net/optim/function/leasqr.html
%   using references by
%   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.
%   Sam Roweis       http://www.cs.toronto.edu/~roweis/notes/lm.pdf
%   Manolis Lourakis http://www.ics.forth.gr/~lourakis/levmar/levmar.pdf
%   Hans Nielson     http://www2.imm.dtu.dk/~hbn/publ/TR9905.ps
%   Mathworks        optimization toolbox reference manual
%   K. Madsen, H.B., Nielsen, and O. Tingleff
%   http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf

global lm_iteration lm_func_calls
% iteration counter
lm_iteration  = 0;
% running count of function evaluations
lm_func_calls = 0;

if nargin == 0
    swhelp ndbase.lm
    return
end

% number of parameters
Np    = numel(p0);
% vertical parameter vector
p0    = p0(:);

inpForm.fname  = {'dp'   'lb'       'ub'      'MaxIter' 'eps1' 'TolX' 'MaxFunEvals' 'confLev'      'win'     };
inpForm.defval = {1e-3   -inf(1,Np) inf(1,Np) 10*Np     1e-3   1e-3   100*Np        erf(1/sqrt(2)) [-inf inf]};
inpForm.size   = {[1 -1] [1 Np]     [1 Np]    [1 1]     [1 1]  [1 1]  [1 1]         [1 1]          [1 2]     };

inpForm.fname  = [inpForm.fname  {'eps2' 'eps3' 'lambda0' 'lUp' 'lDown' 'update' 'extraStat' 'vary'      'nu0'}];
inpForm.defval = [inpForm.defval {1e-2    0.1    1e-2      30    7       'nielsen' true       true(1,Np) 10   }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]  [1 1]     [1 1] [1 1]   [1 -2]    [1 1]      [1 Np]      [1 1]}];

param = sw_readparam(inpForm, varargin{:});

switch param.update
    case 'lm'
        update = 1;
    case 'quadratic'
        update = 2;
    case 'nielsen'
        update = 3;
    otherwise
        error('lm:WrongInput','Wrong ''update'' option!')
end

% check input function
if ischar(func0)
    % convert to fuction handle
    func0 = str2func(func0);
end

if ~isa(func0,'function_handle')
    error('lm:WrongInput','The input function is neither a string, not function handle!');
end

if isempty(dat)
    % define fake data if is not given
    dat   = struct; 
    dat.x = (1:Np)';
    dat.y = dat.x*0;
    dat.e = dat.x*0+1;
    
    func = @(x,p)func0(p)*ones(Np,1);
else
    % just use the given function
    func = func0;
    % make column vectors of all input
    dat.x = dat.x(:);
    dat.y = dat.y(:);
    
    % remove unnecessary data
    xKeep = dat.x>param.win(1) & dat.x<param.win(2);
    dat.x = dat.x(xKeep);
    dat.y = dat.y(xKeep);
    if isfield(dat,'e')
        dat.e = dat.e(xKeep);
    end
end

% TODO: lb = -100*abs(p0); ub = 100*abs(p0)

% number of independent variables
Nx    = numel(dat.x);

% define fit weights
if ~isfield(dat,'e') || isempty(dat.e) || ~any(dat.e(:))
    weight = 1./abs(dat.y);
else
    if any(dat.e(:)<0)
        error('lm:WrongInput','Standard deviations have to be positive!')
    end
    weight = 1./dat.e(:).^2;
end

% previous set of parameters
pOld  = p0(:)';
% previous model values: yOld = func(x,pOld)
yOld  = func(dat.x,pOld);
% empty Jacobian matrix
J     = zeros(Nx,Np);
% empty output
stat = struct('msg',[],'warning',[],'p',[],'sigP',[],'redX2',[],'Rsq',[],'sigY',[],'corrP',[]);

if numel(dat.x) ~= numel(dat.y)
    error('lm:WrongInput',['The number of independent vriables (x) is '...
        'not equal to the number of dependent variables (y)!']);
end

% column vectors of parameters
param.lb = param.lb(:);
param.ub = param.ub(:);

if numel(param.dp) == 1
    param.dp = repmat(param.dp,[1 Np]);
end

if numel(param.vary) == 1
    param.vary = repmat(param.vary,[1 Np]);
end

% fixed the requested parameters
param.vary = logical(param.vary);
if any(~param.vary)
    param.dp(~param.vary) = 0;
end

% indices of the parameters to be fit
idxP = find(param.dp ~= 0);
% number of parameters to vary
Nv = numel(idxP);
% statistical degrees of freedom
DoF   = Nx - Nv + 1;

% termination flag
isFinished = 0;

% initialize Jacobian with finite difference calculation
[JtWJ,JtWdy,X2,yHat,J] = linfitmat(func,dat.x,pOld,yOld,1,J,p0,dat.y,weight,param.dp);

if all(param.dp==0)
    exitFlag   = 7;
    isFinished = 1;
end

if max(abs(JtWdy)) < param.eps1
    exitFlag   = 1;
    isFinished = 1;
end

switch update
    case 1
        % Marquardt: initial lambda
        lambda  = param.lambda0;
    otherwise
        % Quadratic and Nielsen
        lambda  = param.lambda0 * max(diag(JtWJ));
        nu = param.nu0;
end

% default dy value
dy = dat.y*0;

% previous value of X2
X2old = X2;
if param.extraStat
    % initialize convergence history
    cvgHst = ones(0,Np+3);
end

while ~isFinished && lm_iteration <= param.MaxIter
    
    lm_iteration = lm_iteration + 1;
    
    % incremental change in parameters
    switch update
        case 1
            % Marquardt
            h = ( JtWJ + lambda*diag(diag(JtWJ)) ) \ JtWdy;
        otherwise
            % Quadratic and Nielsen
            h = ( JtWJ + lambda*eye(Nv) ) \ JtWdy;
    end
    
    % this is a big step
    %  big = max(abs(h./p)) > 2;
    
    % are parameters [p+h] much better than [p] ?
    % update the [idx] elements
    pTry = p0;
    pTry(idxP) = pTry(idxP) + h;
    % apply constraints
    pTry = min(max(param.lb,pTry),param.ub);
    % residual error using p_try
    dy = dat.y - func(dat.x,pTry');
    % floating point error, break
    if ~all(isfinite(dy))
        exitFlag = 6;
        break
    end
    
    lm_func_calls = lm_func_calls + 1;
    % Chi-squared error criteria
    X2try = dy' * (dy .* weight);
    
    if update == 2
        % Quadratic
        % One step of quadratic line update in the h direction for minimum X2
        alpha =  JtWdy'*h / ( (X2try - X2)/2 + 2*JtWdy'*h ) ;
        h = alpha * h;
        % update only [idx] elements
        pTry = p0;
        pTry(idxP) = pTry(idxP) + h;
        % apply constraints
        pTry = min(max(param.lb,pTry),param.ub);
        % residual error using p_try
        dy = dat.y - func(dat.x,pTry');
        lm_func_calls = lm_func_calls + 1;
        X2try = dy' * (dy.*weight);   % Chi-squared error criteria
    end
    
    rho = (X2 - X2try) / ( h' * (lambda * h + JtWdy) );
    
    if rho > param.eps3
        % it IS significantly better
        dX2   = X2 - X2old;
        X2old = X2;
        pOld  = p0;
        yOld  = yHat;
        % accept p_try
        p0    = pTry(:);
        
        [JtWJ,JtWdy,X2,yHat,J] = linfitmat(func,dat.x,pOld,yOld,dX2,J,p0,dat.y,weight,param.dp);
        
        % decrease lambda ==> Gauss-Newton method
        switch update
            case 1
                % Levenberg
                lambda = max(lambda/param.lDown,1e-7);
            case 2
                % Quadratic
                lambda = max( lambda/(1 + alpha) , 1e-7 );
            case 3
                % Nielsen
                lambda = lambda*max(1/3,1-(2*rho-1)^3);
                nu = param.nu0;
        end
    else
        % it IS NOT better
        % do not accept p_try
        X2 = X2old;
        
        if ~rem(lm_iteration,2*Np)
            % rank-1 update of Jacobian
            [JtWJ,JtWdy,~,yHat,J] = linfitmat(func,dat.x,pOld,yOld,-1,J,p0,dat.y,weight,param.dp);
        end
        
        % increase lambda  ==> gradient descent method
        switch update
            case 1
                % Levenberg
                lambda = min(lambda*param.lUp,1e7);
            case 2
                % Quadratic
                lambda = lambda + abs((X2try - X2)/2/alpha);
            case 3
                % Nielsen
                lambda = lambda * nu;
                nu     = nu*param.nu0;
        end
        
    end
    
    if param.extraStat
        % update convergence history ... save _reduced_ Chi-square
        cvgHst(lm_iteration,:) = [lm_func_calls  p0'  X2/DoF lambda];
    end
    
    if max(abs(JtWdy)) < param.eps1 && lm_iteration > 2
        exitFlag = 1;
        isFinished = true;
    end
    if max(abs(h./p0(idxP))) < param.TolX  &&  lm_iteration > 2
        exitFlag = 2;
        isFinished = true;
    end
    if X2/DoF < param.eps2 &&  lm_iteration > 2
        exitFlag = 3;
        isFinished = true;
    end
    if lm_iteration >= param.MaxIter
        exitFlag = 4;
        isFinished = true;
    end
    
    if lm_func_calls >= param.MaxFunEvals
        exitFlag = 5;
        isFinished = true;
    end
    
end

switch exitFlag
    case 1
        stat.msg = 'Convergence in r.h.s. ("JtWdy").';
    case 2
        stat.msg = 'Convergence in parameters (dX<TolX).';
    case 3
        stat.msg = 'Convergence in reduced Chi-square.';
    case 4
        stat.msg = 'Maximum Number of iterations is reached without convergence, increase ''MaxIter''!';
    case 5
        stat.msg = 'Maximum Number of function evaluations is reached without convergence, increase ''MaxFunEvals''!';
    case 6
        stat.msg = 'Floating point error, parameter change is smaller than eps!';
    case 7
        stat.msg = 'All parameters are fixed!';
end

if exitFlag > 3 && exitFlag < 7
    stat.error = true;
    warning('lm:convergence',stat.msg)
else
    stat.warning = false;
end

% convergence achieved, find covariance and confidence intervals
if var(weight) == 0
    % recompute equal weights for paramter error analysis
    weight = DoF/(dy'*dy) * ones(Nx,1);
end

if nargout > 1
    % reduced Chi-square
    stat.redX2 = X2/DoF;
end

[JtWJ,~,~,yHat,J0] = linfitmat(func,dat.x,pOld,yOld,-1,J,p0,dat.y,weight,param.dp);

% extend matrix to original parameter size
J = zeros(Nx,Np);
J(:,idxP) = J0;

% standard error of parameters
covP0 = inv(JtWJ);
covP = zeros(Np,Np);
covP(idxP,idxP) = covP0;

stat.sigP = sqrt(diag(covP))';

if param.extraStat
    
    % error of the fit at the given confidence level
    stat.sigY = sqrt(2)*erfinv(param.confLev)*sqrt(diag(J*covP*J'))';
    
    % parameter correlation matrix
    stat.corrP = covP./(stat.sigP*stat.sigP');
    
    % coefficient of multiple determination
    %stat.Rsq = corr([dat.y yHat]);
    stdy = std(dat.y)*std(yHat);
    if stdy == 0
        stdy = 1;
    end
    
    stat.Rsq = cov(dat.y,yHat)./stdy;
    stat.Rsq = stat.Rsq(1,2).^2;
    
    % convergence history
    stat.cvgHst = cvgHst;
    
end

% parameter values
stat.p = p0';
% save counters
stat.nIter      = lm_iteration;
stat.nFunEvals  = lm_func_calls;
stat.algorithm  = 'Levenberg-Marquardt';
stat.func       = func;
stat.exitFlag   = exitFlag;

% save the input parameters
stat.param = param;

% type of update
switch update
    case 1
        stat.param.updateStr = 'Levenberg-Marquardt lambda update';
    case 2
        stat.param.updateStr = 'Quadratic update';
    case 3
        stat.param.updateStr = 'Nielsen''s lambda update equations';
end

% final model value
yHat = yHat';

% final parameter value
p0 = p0';

end


function J = jacobian(func,x,p,y,dp)
% numerical partial derivative calculator
%
% J = JACOBIAN(func,x,p,y,dp)
%
% The function calculates the partial derivatives (Jacobian) dfunc(x,p)/dp
% via finite differences. Requires n or 2n function evaluations, depending
% on the method.
%
% Input:
% func  Function handle with the following header:
%           y = func(x,p)
%       where y is a vector and we derivate according to every element of
%       p.
% x     Vector of N independent variables.
% p  	Vector of M parameters.
% y     Vector of N precalculated function values.
% dp 	Fractional increment of p for numerical derivatives:
%           dp(j)>0     central differences calculated
%           dp(j)<0     one sided differences calculated
%           dp(j)=0     sets corresponding partials to zero, i.e. holds p(j)
%                       fixed.
%       Default value is 1e-3.
%
% Output:
%
% J  	Jacobian matrix with dimensions [N,M] and values:
%           J(i,j) = dFunc(x,p)(i) / dp(j)
%

%   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. November 2005
%   modified from: ftp://fly.cnuce.cnr.it/pub/software/octave/leasqr/
%   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.

global  lm_func_calls

% number of data points
M = numel(y);
% indices of parameters to vary
pIdx = find(dp~=0);

% number of varied parameters
Nv = numel(pIdx);

ps  = p;
J   = zeros(M,Nv);
% initialize Jacobian to Zero
del = zeros(Nv,1);

% loop over all parameters
for ii = 1:Nv
    idx = pIdx(ii);
    % parameter perturbation
    del(idx) = dp(idx) * (1+abs(p(idx)));
    % perturb parameter p(j)
    p(idx)   = ps(idx) + del(idx);
    
    if del(idx)%~=0
        y1 = func(x,p');
        lm_func_calls = lm_func_calls + 1;
        
        if dp(idx) < 0
            % backwards difference
            J(:,ii) = (y1-y)/del(idx);
        else
            % central difference, additional func call
            p(idx)   = ps(idx) - del(idx);
            J(:,ii) = (y1-func(x,p'))/(2*del(idx));
            lm_func_calls = lm_func_calls + 1;
        end
    end
    % restore p(j)
    p(idx) = ps(idx);
    
end

end

function J = broyden(pOld,yOld,J,p,y,dp)
% carry out a rank-1 update to the Jacobian matrix using Broyden's equation
%
% J = BROYDEN(p_old,y_old,J,p,y,dp)
%
% Input:
%
% pOld 	Previous set of parameters.
% yOld 	Model evaluation at previous set of parameters, yHat(x,pOld)
% J  	Current value of the Jacobian matrix.
% p     Current  set of parameters.
% y     Model evaluation at current  set of parameters, yHat(x,p)
% dp    Only determine the jacobian where dp~=0.
%
% Output:
%
% J  	Jacobian matrix with dimensions [N,M] and values:
%           J(i,j) = dFunc(x,p)(i) / dp(j)
%

h = p(dp~=0) - pOld(dp~=0);
% Broyden rank-1 update eq'n
J = J + (y-yOld-J*h)*h'/(h'*h);

end


function [JtWJ,JtWdy,chiSq,yFunc,J] = linfitmat(func,x,pOld,yOld,dX2,J,p,yDat,weight,dp)
% evaluate the linearized fitting matrix
%
% [JtWJ,JtWdy,chiSq,yHat,J] = LINFITMAT(func,x,pOld,yOld,dX2,J,p,yDat,weight,dp)
%
% Evaluates the linearized fitting matrix, JtWJ, and vector JtWdy, and
% calculate the Chi-squared error function Chi_sq.
%
% Input:
%
% func  Function handle with the following header:
%           y = func(x,p)
%       where y is a vector and we derivate according to every element of
%       p.
% x     Vector of N independent variables.
% pOld  Vector of M parameters of previous values.
% yOld  Vector of N previous calculated model values.
% dX2   Previous change in Chi-squared criteria.
% J  	Jacobian matrix with dimensions [N,M] and values:
%           J(i,j) = dFunc(x,p)(i) / dp(j)
% p  	Vector of M parameters.
% y     Vector of N precalculated function values.
% yDat  Vector of data to be fitted.
% weight The weighting vector for least squares fit ...
% dp 	Fractional increment of p for numerical derivatives:
%           dp(j)>0     central differences calculated
%           dp(j)<0     one sided differences calculated
%           dp(j)=0     sets corresponding partials to zero, i.e. holds p(j)
%                       fixed.
%       Default value is 1e-3.
%
% Output:
%
% JtWJ      Linearized Hessian matrix (inverse of covariance matrix).
% JtWdy   	Linearized fitting vector.
% chiSq     Chi-squared criteria, weighted sum of the squared residuals.
% yHat  	Model evaluated with parameters p.
% J         Jacobian matrix.
%

%   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. November 2005
%   modified from: ftp://fly.cnuce.cnr.it/pub/software/octave/leasqr/
%   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.

global   lm_iteration  lm_func_calls

% number of varying parameters
Np = numel(p);
% number of varying parameters
Nv = sum(dp~=0);

% evaluate model using parameters 'p'
yFunc = func(x,p');
lm_func_calls = lm_func_calls + 1;

if ~rem(lm_iteration,2*Np) || dX2 > 0
    % finite difference
    J = jacobian(func,x,p,yFunc,dp);
else
    % rank-1 update
    J = broyden(pOld,yOld,J,p,yFunc,dp);
end

% residual error between model and data
dy = yDat - yFunc;
% Chi-squared error criteria
chiSq = dy' * (dy.* weight);

JtWJ  = J' * (J.*repmat(weight,[1,Nv]));

JtWdy = J' * (weight.*dy);

end