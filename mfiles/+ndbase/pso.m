function [pOpt,fVal,stat] = pso(dat,func,p0,varargin)
% minimization via particle swarm optimisation
%
% [pOpt,fVal,stat] = NDBASE.PSO([],func,p0,'Option1','Value1',...)
%
% [pOpt,fVal,stat] = NDBASE.PSO(dat,func,p0,'Option1','Value1',...)
%
%
% PSO finds a minimum of a function of several variables using the particle
% swarm optimization (PSO) algorithm originally introduced in 1995 by
% Kennedy and Eberhart. This algorithm was extended by Shi and Eberhart in
% 1998 through the introduction of inertia factors to dampen the velocities
% of the particles. In 2002, Clerc and Kennedy introduced a constriction
% factor in PSO, which was later on shown to be superior to the inertia
% factors. Therefore, the algorithm using a constriction factor was
% implemented here.
%
% The function has two different modes, depending on the first input
% argument. If dat is empty, pso minimizes the cost function func, that has
% the following header:
%                           R2 = func(p)
%
% If dat is a struct, the function optimizes the model defined by func via
% least squares to the data stored in dat. In this case func has the
% following header:
%                           yModel = func(x,p)
% And the least square deviation is defined by:
%                       R2 = sum((yData-yModel).^2/errorData.^2)
%
%
% Input:
%
% dat       Either empty or contains data to be fitted stored in a
%           structure with fields:
%               dat.x   vector of N independent variables,
%               dat.y   vector of N data values to be fitted,
%               dat.e   vector of N standard deviation (positive numbers)
%                       used to weight the fit. If zero or missing
%                       1./dat.y^2 will be assigned to each point.
% func      Function handle with one of the following definition:
%               R2 = func(p)        if dat is empty,
%               y  = func(x,p)      if dat is a struct.
%           Here x is a vector of N independent variables, p are the
%           M parameters to be optimized and y is the simulated model, R2
%           is the value to minimize.
% p0        Vector of M initial parameters.
%
% Options:
%
% lb        Vector with N elements, lower boundary of the parameters.
%           Default value is -1e5.
% ub        Vector with N elements, upper boundary of the parameters.
%           Default value is 1e5.
% MaxIter   Maximum number of iterations, default value is 100*M.
% MaxFunEvals Maximum number of function evaluations, default value is
%           1000*M.
% TolX      Convergence tolerance for parameters, default value is 1e-3.
% Display   Level of information to print onto the Command Line during
%           optimisation. Possible values are 'off', 'notify' and 'iter'.
%           Default is 'off'.
% TolFun    Convergence tolerance on the R2 value (return value of func or
%           the weighted least square deviation from data). Default value
%           is 1e-3.
% SwarmC1   Swarm parameter, also called phi_p in the literature. Default
%           value is 2.8.
% SwarmC2   Swarm parameter, also called phi_g in the literature. Default
%           value is 1.3.
% k0        Swarm paramter, it can take values between 0 and 1, but is
%           usually set to one (Montes de Oca et al., 2006). Default value
%           is 1.
% seed      Seed number for the random number generator, by default it is
%           seeded from the system clock every time the function is called.
% PopulationSize Size of the swarm, default value is 25.
% autoTune  Determine the size of the swarm based on the number of free
%           parameters (number of parameters that are not fixed by setting
%           lb(i)==ub(i)). Default is false. The optimal swarm size is:
%               PopulationSize = 25 + 1.4*Nf
%           This value is based on the technical report:
%               Good Parameters for Particle Swarm Optimization By Magnus Erik Hvass Pedersen Hvass Laboratories
%               Technical Report no. HL1001
%               http://www.hvass-labs.org/people/magnus/publications/pedersen10good-pso.pdf
%
% See also NDBASE.LM.
%

% Copyright (C) 2006 Brecht Donckels, BIOMATH, brecht.donckels@ugent.be
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

if nargin == 0
    help ndbase.pso
    return
end

% number of parameters
Np    = numel(p0);

oNp = ones(1,Np);

inpForm.fname  = {'Display' 'TolFun' 'TolX' 'MaxIter' 'MaxFunEvals' 'SwarmC1' 'seed'       };
inpForm.defval = {'off'     1e-3     1e-3   100*Np    1000*Np       2.8      sum(100*clock)};
inpForm.size   = {[1 -1]    [1 1]    [1 1]  [1 1]     [1 1]         [1 1]    [1 1]         };

inpForm.fname  = [inpForm.fname  {'SwarmC2' 'PopulationSize' 'lb'      'ub'     'autoTune' 'k0'}];
inpForm.defval = [inpForm.defval {1.3       25               -1e5*oNp   1e5*oNp false      1   }];
inpForm.size   = [inpForm.size   {[1 1]     [1 1]            [1 Np]     [1 Np]  [1 1]     [1 1]}];

param = sw_readparam(inpForm, varargin{:});

% parameter boundaries
UB = param.ub;
LB = param.lb;

if any(UB<LB)
    error('pso:WrongInput','Upper boundary has to be larger than the lower boundary!');
end

% number of free parameters
Nf = sum(UB>LB);

if param.autoTune
    param.PopulationSize = 25 + 1.4*Nf;
end

% population size
popSize = param.PopulationSize;

% maximum number of iterations
maxIter = param.MaxIter;

% swarm parameters
swarmC1 = param.SwarmC1;
swarmC2 = param.SwarmC2;

% check input function
if ischar(func)
    % convert to fuction handle
    func = str2func(func);
end

if ~isa(func,'function_handle')
    error('pso:WrongInput','The input function is neither a string, not function handle!');
end

% define weighted least squares if dat is given
if ~isempty(dat)
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
% TODO

% set EXITFLAG to default value
exitFlag = 0;

% seed the random number generator
rng(param.seed);

% initialize swarm (each row of swarm corresponds to one particle)
swarm = zeros(popSize,Np,maxIter);

swarm(1,:,1) = p0(:)';
for ii = 2:popSize
    swarm(ii,:,1) = LB(:)' + rand(1,Np).*(UB(:)'-LB(:)');
end

% initialize parameters
velocities          = zeros(popSize,Np,maxIter);
fitness             = nan(popSize,maxIter);
pBest               = nan(popSize,Np,maxIter);
gBest               = nan(maxIter,Np);
pBestFitness        = nan(popSize,maxIter);
gBestFitness        = nan(popSize,1);
indexPBest          = nan(popSize,maxIter);
indexGBestParticle  = nan(maxIter,1);
indexGBestIteration = nan(maxIter,1);

% calculate constriction factor from acceleration coefficients
if swarmC1+swarmC2 <= 4
    % Display message
    if strcmp(param.Dispplay,'iter') || strcmp(param.Display,'notify')
        fprintf(['Sum of Cognitive Acceleration Coefficient and Social Acceleration Coefficient is less then or equal to 4.\n'...
            'Their values were adjusted automatically to satisfy this condition.\n\n']);
    end
    % the values are adjusted so that the sum is equal to 4.1, keeping the ratio SwarmC1/SwarmC2 constant
    swarmC1 = swarmC1*4.1/(swarmC1+swarmC2);
    swarmC2 = swarmC2*4.1/(swarmC1+swarmC2);
end

% calculate constriction factor
param.ConstrictionFactor = 2*param.k0/(abs(2-(swarmC1+swarmC2)-sqrt((swarmC1+swarmC2)^2-4*(swarmC1+swarmC2))));

% for each iteration....
for ii = 1:maxIter
    
    % calculate FITNESS values for all particles in SWARM
    % (each row of FITNESS corresponds to the FITNESS value of one particle)
    % (each column of FITNESS corresponds to the FITNESS values of the particles in one iteration)
    
    for jj = 1:popSize
        fitness(jj,ii) = calculate_cost(func,swarm(jj,:,ii),LB,UB);
    end
    
    % identify particle's location at which the best FITNESS has been achieved (PBEST)
    for jj = 1:popSize
        [pBestFitness(jj,ii),indexPBest(jj,ii)] = min(fitness(jj,:));
        pBest(jj,:,ii) = swarm(jj,:,indexPBest(jj,ii));
    end
    
    % identify the particle from the SWARM at which the best FITNESS has been achieved so far (GBEST)
    [gBestFitness(ii),index_gbest] = min(reshape(fitness,numel(fitness),1));
    [indexGBestParticle(ii),indexGBestIteration(ii)] = ind2sub(size(fitness),index_gbest);
    gBest(ii,:) = swarm(indexGBestParticle(ii),:,indexGBestIteration(ii));
    
    % update the VELOCITIES
    velocities(:,:,ii+1) = param.ConstrictionFactor.*(velocities(:,:,ii) + ...
        swarmC1.*rand(popSize,Np).*(pBest(:,:,ii)-swarm(:,:,ii)) + ...
        swarmC2.*rand(popSize,Np).*(repmat(gBest(ii,:),[popSize 1 1 1])-swarm(:,:,ii)));
    
    % update particle positions
    swarm(:,:,ii+1) = swarm(:,:,ii)+velocities(:,:,ii+1);
    
    % to make sure that particles stay within specified bounds...
    %   (suppose that the particle's new position is outside the boundaries,
    %    then the particle's position is adjusted by assuming that the boundary
    %    acts like a wall or mirror) (selfmade solution)
    for jj = 1:popSize
        for kk = 1:Np
            % check upper boundary
            if numel(UB) == 1
                if swarm(jj,kk,ii+1) > UB
                    swarm(jj,kk,ii+1) = UB-rand*(swarm(jj,kk,ii+1)-UB);
                    velocities(jj,kk,ii+1) = swarm(jj,kk,ii+1)-swarm(jj,kk,ii);
                end
            else
                if swarm(jj,kk,ii+1) > UB(kk)
                    swarm(jj,kk,ii+1) = UB(kk)-rand*(swarm(jj,kk,ii+1)-UB(kk));
                    velocities(jj,kk,ii+1) = swarm(jj,kk,ii+1)-swarm(jj,kk,ii);
                end
            end
            % check lower boundary
            if numel(UB) == 1
                if swarm(jj,kk,ii+1) < LB
                    swarm(jj,kk,ii+1) = LB+rand*(LB-swarm(jj,kk,ii+1));
                    velocities(jj,kk,ii+1) = swarm(jj,kk,ii+1)-swarm(jj,kk,ii);
                end
            else
                if swarm(jj,kk,ii+1) < LB(kk)
                    swarm(jj,kk,ii+1) = LB(kk)+rand*(LB(kk)-swarm(jj,kk,ii+1));
                    velocities(jj,kk,ii+1) = swarm(jj,kk,ii+1)-swarm(jj,kk,ii);
                end
            end
        end
    end
    
    % give user feedback on screen if requested
    if strcmp(param.Display,'iter')
        if ii == 1
            fprintf(' Nr Iter  Nr Fun Eval    Current best function    Current worst function       Best function\n');
        end
        fprintf(' %5.0f %5.0f %12.6g %12.6g %15.6g\n',ii,ii*popSize,min(fitness(:,ii)),max(fitness(:,ii)),gBestFitness(ii));
    end
    
    if param.TolX >0 && max(max(abs(diff(swarm(:,:,ii),1,1)))) < param.TolX
        % maximum difference between the coordinates of the vertices in simplex is less than TolX
        stat.msg  = 'Convergence in parameters (dX<TolX).';
        exitFlag = 2;
    end
    if param.TolFun && abs(min(fitness(:,ii))-gBestFitness(ii)) < abs(param.TolFun) ...
            && abs(min(fitness(:,ii))-gBestFitness(ii)) > 0
        % difference between best and worst function evaluation in simplex is smaller than TolFun
        exitFlag = 7;
        stat.msg  = 'Convergence in function value change (dF<TolFun).';
    end
    
    % if a termination criterium was met, the value of EXITFLAG should have changed
    % from its default value of -2 to -1, 0, 1 or 2
    if exitFlag
        break
    end
    
end

if exitFlag == 0
    % no convergence,but maximum number of iterations has been reached
    stat.msg  = 'Maximum Number of iterations is reached without convergence!';
    exitFlag   = 4;
    stat.warning = true;
    warning('pso:convergence','Convergence is not reached!')
else
    stat.warning = false;
end

% return solution
pOpt = gBest(ii,:);
fVal = gBestFitness(ii);

% store number of iterations
stat.p          = pOpt;
stat.sigP       = [];
if isempty(dat)
    stat.redX2 = fVal;
else
    % divide R2 with the statistical degrees of freedom
    stat.redX2   = fVal/(numel(dat.x)-Np+1);
end

stat.Rsq        = [];
stat.sigY       = [];
stat.corrP      = [];
stat.cvgHst     = [];
stat.nIter      = ii;
stat.nFunEvals  = ii*popSize;
stat.algorithm  = 'Particle Swarm Optimization';
if isempty(dat)
    stat.func   = func;
else
    stat.func   = func0;
end

stat.exitFlag   = exitFlag;
stat.param      = param;

% store number of function evaluations
if ~isempty(dat)
    fVal = func0(dat.x,pOpt);
end

end

function [YTRY] = calculate_cost(FUN,PTRY,LB,UB)
% cost function evaluation

for ii = 1:numel(LB)
    % check lower bounds
    if PTRY(ii) < LB(ii)
        YTRY = 1e12+(LB(ii)-PTRY(ii))*1e6;
        return
    end
    % check upper bounds
    if PTRY(ii) > UB(ii)
        YTRY = 1e12+(PTRY(ii)-UB(ii))*1e6;
        return
    end
end

% calculate cost associated with PTRY
YTRY = FUN(PTRY);

end