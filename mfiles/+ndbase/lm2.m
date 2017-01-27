function [pOpt, fVal, stat] = lm2(dat,func0,p0,varargin)
% optimization of parameters using the Levenberg Marquardt method
%
% [pOpt,fVal,stat] = NDBASE.LM([],func,p0,'Option1','Value1',...)
%
% Levenberg Marquardt curve-fitting function, minimizes the return value of
% a given function.
%
% Input:
%
% func      Function handle with the following definition:
%               R = func(p)
%           where p are the M parameters to be optimized.
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

% Reference:
% Fletcher, R., (1971): A Modified Marquardt Subroutine for Nonlinear Least
% Squares. Rpt. AERE-R 6799, Harwell

% Miroslav Balda,
% balda AT cdm DOT cas DOT cz
% 2007-07-02    v 1.0
% 2008-12-22    v 1.1 * Changed name of the function in LMFsolv
%                     * Removed part with wrong code for use of analytical
%                       form for assembling of Jacobian matrix
% 2009-01-08    v 1.2 * Changed subfunction printit.m for better one, and
%                       modified its calling from inside LMFsolve.
%                     * Repaired a bug, which caused an inclination to
%                       istability, in charge of slower convergence.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Copyright (c) 2007, Miroslav Balda
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are
%  met:
%
%      * Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in
%        the documentation and/or other materials provided with the distribution
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.

p  = p0(:)';
Np = numel(p);

inpForm.fname  = {'dp'   'lb'       'ub'      'MaxIter' 'scaleD' 'TolX' 'MaxFunEvals' 'confLev'      'win'     };
inpForm.defval = {1e-3   -inf(1,Np) inf(1,Np) 100*Np    'auto'   1e-4   1000*Np       erf(1/sqrt(2)) [-inf inf]};
inpForm.size   = {[1 -1] [1 Np]     [1 Np]    [1 1]     [1 -2]    [1 1]  [1 1]         [1 1]          [1 2]     };


inpForm.fname  = [inpForm.fname  {'TolFun' 'eps3' 'lambda0' 'vary'      'nu0'}];
inpForm.defval = [inpForm.defval {1e-6     0.1    1e-2      true(1,Np)  10   }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]  [1 1]     [1 Np]      [1 1]}];

param = sw_readparam(inpForm, varargin{:});

% check input function
if ischar(func0)
    % convert to fuction handle
    func0 = str2func(func0);
end

if ~isa(func0,'function_handle')
    error('lm:WrongInput','The input function is neither a string, not function handle!');
end

% define weighted least squares if dat is given
if ~isempty(dat)
    dat.x = dat.x(:);
    dat.y = dat.y(:);
    
    if ~isfield(dat,'e') || isempty(dat.e) || ~any(dat.e(:))
        func = @(p)sum((func0(dat.x(:),p)-dat.y(:)).^2);
    else
        if any(dat.e(:)<0)
            error('pso:WrongInput','Standard deviations have to be positive!')
        end
        weight = 1./dat.e(:).^2;
        func = @(p)sum(weight.*(func0(dat.x(:),p)-dat.y(:)).^2);
    end
    
else
    func = func0;
end

% residuals at starting point
resid = sqrt(func(p));
nFunEvals = 1;
if numel(resid) == 1
    resid = resid*ones(5,1)/5;
else
    resid = resid(:);
end

S    = resid'*resid;
epsx = param.TolX(:);
J    = finjac(func,resid,p,epsx);
nFunEvals = nFunEvals + Np;

% system matrix
A = J.'*J;
v = J.'*resid;

if strcmp(param.scaleD,'auto')
    % automatic scaling
    D = diag(A);
    D(D(1:Np)==0)=1;
    D=diag(D);
else
    D = param.scaleD;
    if numel(D)>1
        D = diag(sqrt(abs(D(1:Np)))); % vector of individual scaling
    else
        D = sqrt(abs(D))*eye(Np);     % scalar of unique scaling
    end
end

Rlo = 0.25;
Rhi = 0.75;
l   = 1;
lc  = 0.75;
cnt = 0;
exitflag = 0;

while exitflag == 0      %   MAIN ITERATION CYCLE
    try
        d = pinv(A+l*D)*v;
    catch
        % negative solution increment
        d = (A+l*D)\v;
    end
    newP = p-d';
    rd   = sqrt(func(newP));
    nFunEvals = nFunEvals+1;
    
    if length(rd) == 1
        rd=rd*ones(5,1)/5;
    else
        rd=rd(:);
    end
    
    
    Sd = rd.'*rd;
    % predicted reduction
    dS = d.'*(2*v-A*d);
    R=(S-Sd)/dS;
    
    if R>Rhi
        l=l/2;
        if l<lc
            l=0;
        end
    elseif R<Rlo
        nu=(Sd-S)/(d.'*v)+2;
        if nu<2
            nu=2;
        elseif nu>10
            nu=10;
        end
        if l==0
            lc=1/max(abs(diag(pinv(A))));
            l=lc;
            nu=nu/2;
        end
        l=nu*l;
    end
    cnt=cnt+1;
    
    S   = Sd;
    p   = newP;
    resid   = rd;
    J   = finjac(func,resid,p,epsx);
    nFunEvals = nFunEvals+Np;
    
    A = J.'*J;
    v = J.'*resid;
    
    
    if cnt >= param.MaxIter
        % max iteration reached
        exitflag = -2;
    elseif all(abs(d)<param.TolX/1000)
        % parameter change increment is negligible
        exitflag = -1;
    elseif all(abs(resid)<param.TolFun)
        % function change increment reached
        exitflag = -5;
    elseif nFunEvals > param.MaxFunEvals
        % max nb function evaluations reached
        exitflag = -3;
    end
end

%   final solution
pOpt = p;

stat.S         = S;
stat.exitflag  = exitflag;
fVal           = func(pOpt);
stat.nFunEvals = nFunEvals;

end

function J = finjac(func,r,p,epsx)
% numerical approximation to Jacobi matrix
%

% pars=column, function=row vector or scalar
Np = numel(p);
J  = zeros(numel(r), Np);
p  = p(:)';
r  = r(:);
if numel(epsx)==1
    epsx=epsx*max(abs(p),1);
end
if any(epsx == 0)
    epsx(~epsx) = 1e-4;
end

epsx=epsx(:)';

for kk = 1:Np
    dx = 0.25*epsx(kk);
    xd = p;
    xd(kk) = xd(kk)+dx;
    rd = sqrt(func(xd));
    if numel(rd) == 1
        rd=rd*ones(5,1)/5;
    else
        rd=rd(:);
    end
    
    if dx>0
        J(:,kk)=((rd-r)/dx);
    end
end

end