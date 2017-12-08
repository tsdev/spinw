function fitsp = fitpow(obj, data, Q_slices, varargin)
% fits spin wave spectra to experimental spectral data
%
% fitsp = FITPOW(obj, 'Option1', Value1, ...)
%
% Options:
%
% func      Function to change the Hamiltonian in obj, it has the following
%           header:
%                    obj = @func(obj, x);
% datapath  Path to the file that stores the experimental data. For the
%           input data format see <a href="matlab:doc sw_readspec">sw_readspec</a>.
% Evect     Vector, defines the energy binning of the calculated
%           dispersion. Larger binning steps solve the issue of fitting
%           unresolved modes. Size is [1 nE].
% xmin      Minimum limit of the optimisation parameters, optional.
% xmax      Maximum limit of the optimisation parameters, optional.
% x0        Starting value of the optimisation parameters. If empty
%           or undefined, then random values are used.
% nRun      Number of consecutive fitting runs, each result is saved in the
%           output fitsp.x and R arrays. If the Hamiltonian given by the
%           random x parameters is incompatible with the ground state,
%           those x values will be skipped and new random x values will be
%           generated. Default is 1.
% nMax      Maximum number of runs, including the ones that produce error
%           (due to incompatible ground state). Default is 1000.
% hermit    Method for matrix diagonalization:
%                  true      J.H.P. Colpa, Physica 93A (1978) 327,
%                  false     R.M. White, PR 139 (1965) A450.
%           Colpa: the grand dynamical matrix is converted into another
%                  Hermitian matrix, that will give the real eigenvalues.
%           White: the non-Hermitian g*H matrix will be diagonalised,
%                  that is not the elegant method.
%           Advise:
%           Always use Colpa's method that is faster, except when small
%           imaginary eigenvalues are expected. In this case only White's
%           method work.
%           Default is true.
% epsilon   Small number that controls wether the magnetic structure is
%           incommensurate or commensurate, default value is 1e-5.
%
% Parameters for visualizing the fit results:
%
% plot      If true, the measured dispersion is plotted together with the
%           fit. Default is true.
% iFact     Factor of the plotted simulated spin wave intensity (red
%           ellipsoids).
% lShift   Vertical shift of the Q point labels on the plot.
%
%
% Optimisation options:
%
% tolx          Minimum change of x when convergence reached, default
%               value is 1e-4.
% tolfun        Minimum change of the R value when convergence reached,
%               default value is 1e-5.
% maxfunevals   Maximum number of function evaluations, default value is
%               1e7.
%
% Output:
%
% Output fitsp is struct type with the following fields:
% obj       Copy of the input sw class object, with the best fitted
%           Hamiltonian.
% x         Final values of the fitted parameters, dimensions are
%           [nRun nPar]. The rows of x are sorted according to increasing R
%           values.
% R         R-value, goodness of the fit, dimensions are [nRun 1], sorted
%           in increasing order.
% exitflag  Exit flag of the fminsearch command.
% output    Output of the fminsearch command.
%
% Any option used by SW.SPINWAVE function are also accepted.
%
% See also SPINW.SPINWAVE, SW_EGRID, SW_NEUTRON, SW_READSPEC, FMINSEARCH.
%

pref = swpref;
tid0 = pref.tid;
T0 = obj.single_ion.T;
x0 = squeeze(mat2cell(tri.matrix.mat,3,3,ones(1,size(tri.matrix.mat,3))))';
is_diag = cellfun(@(x) all(diff(diag(x))==0),x0);
x0(is_diag) = cellfun(@(x) x(1,1),x0(is_diag),'UniformOutput',false);

% These are the powspec options.
inpForm.fname  = {'nRand' 'T'   'formfact' 'formfactfun' 'tid' 'nInt'  'dE'    'x0'      'x0_lab'};
inpForm.defval = {100     T0    false      @sw_mff       tid0  1e3     2       []        []};
inpForm.size   = {[1 1]   [1 1] [1 -2]     [1 1]         [1 1] [1 1]   [1, 1]  [1, -6]   [1, -6]};

inpForm.fname  = [inpForm.fname  {'hermit' 'gtensor' 'title' 'specfun' 'imagChk' }];
inpForm.defval = [inpForm.defval {true     false     title0  @spinwave  true     }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]     [1 -3]  [1 1]      [1 1]    }];

inpForm.fname  = [inpForm.fname  {'extrap' 'fibo' 'optmem' 'binType' 'component'}];
inpForm.defval = [inpForm.defval {false    true   0        'cbin'    'Sperp'    }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]  [1 1]    [1 -4]     [1 -5]    }];

inpForm.fname  = [inpForm.fname  {'fid'}];
inpForm.defval = [inpForm.defval {0   }];
inpForm.size   = [inpForm.size   {[1 1]}];

warning('off','sw_readparam:UnreadInput')
param_PS  = sw_readparam(inpForm, varargin{:});
warning('on','sw_readparam:UnreadInput')

% PSO options.
inpForm.fname  = {'Display' 'TolFun' 'TolX' 'MaxIter' 'SwarmC1' 'seed'        'nRun'};
inpForm.defval = {'off'     1e-3     1e-3   100       2.8      sum(100*clock) 1     };
inpForm.size   = {[1 -1]    [1 1]    [1 1]  [1 1]     [1 1]    [1 1]          [1, 1]};

inpForm.fname  = [inpForm.fname  {'SwarmC2' 'PopulationSize' 'xmin'     'xmax'  'autoTune' 'k0'}];
inpForm.defval = [inpForm.defval {1.3       25               []         []      false      1   }];
inpForm.size   = [inpForm.size   {[1 1]     [1 1]            [1 -2]     [1 -3]  [1 1]      [1 1]}];

warning('off','sw_readparam:UnreadInput')
param_PSO  = sw_readparam(inpForm, varargin{:});
warning('on','sw_readparam:UnreadInput')

% number of parameters (length of x)
nPar = max(max(length(param_PSO.xmin),length(param_PSO.xmax)),length(param_PS.x0));
if (param_PSO.MaxIter)==100
    param_PSO.MaxIter = nPar*param_PSO.MaxIter;
end
nRun = param_PSO.nRun;

% setup parameters
param_PS.tid  = 0;

x     = zeros(nRun,nPar);
redX2 = zeros(nRun,1);

sw_timeit(0,1,param.tid,'Fitting spin wave powder spectra');

idx = 1;
idxAll = 1;

while idx <= nRun
    %try
    if ~isempty(param_PS.x0)
        x0 = param_PS.x0;
    elseif isempty(param_PS.xmax) && isempty(param_PS.xmin)
        x0 = rand(1,nPar);
    elseif isempty(param_PS.xmax)
        x0 = param_PS.xmin + rand(1,nPar);
    elseif isempty(param.xmin)
        x0 = param_PS.x0 - rand(1,nPar);
    else
        x0 = rand(1,nPar).*(param.xmax-param.xmin)+param.xmin;
    end
    
    switch param.optimizer
        case 'pso'
            [x(idx,:),~, output(idx)] = ndbase.pso(dat,@(x,p)spec_fitfun(obj, data, param.func, p, param0),x0,'lb',param.xmin,'ub',param.xmax,...
                'TolX',param.tolx,'TolFun',param.tolfun,'MaxIter',param.maxiter);
            
            redX2(idx) = output.redX2;
            
        case 'simplex'
            [x(idx,:),~, output(idx)] = ndbase.simplex(dat,@(x,p)spec_fitfun(obj, data, param.func, p, param0),x0,'lb',param.xmin,'ub',param.xmax,...
                'TolX',param.tolx,'TolFun',param.tolfun,'MaxIter',param.maxiter);
            
            redX2(idx) = output.redX2;
            
        case 'lm'
            % does not work due to the binning of the spectrum
            % newspec_fitfun is required for this
            [x(idx,:),~, output(idx)] = ndbase.lm(dat,@(x,p)spec_fitfun(obj, data, param.func, p, param0),x0,'lb',param.xmin,'ub',param.xmax,...
                'TolX',param.tolx,'TolFun',param.tolfun,'MaxFunEvals',param.maxfunevals);
            redX2(idx) = output(idx).redX2;
        otherwise
            error('spinw:fitspec:WrongOption','The given optimizer is not supported!')
            
    end
    
    idx = idx + 1;
    %     catch
    %         warning('Hamiltonian is not compatible with the magnetic ground state!');
    %         if ~isempty(param.x0)
    %             error('spinw:fitspec:WrongInput','x0 input incompatible with ground state!');
    %         end
    %
    %         if idxAll >= param.nMax
    %             warning('spinw:fitspec:WrongInput','Maximum number of wrong runs (param.nMax) reached, calculation is terminated!');
    %             break;
    %         end
    %     end
    idxAll = idxAll + 1;
    sw_timeit(idx/nRun*100,0,param.tid);
end

sw_timeit(100,2,param.tid);

% Sort results
[redX2, sortIdx] = sort(redX2);
x = x(sortIdx,:);

% Draw plot of the final result if requested
if param.plot
    figure;
    param0.plot = true;
    % plot the best result if requested
    spec_fitfun(obj, data, param.func, x(1,:), param0);
end

% set the best fit to the spinw object
param.func(obj,x(1,:));

% Store all output in a struct variable.
fitsp.obj      = copy(obj);
fitsp.x        = x;
fitsp.redX2    = redX2;
%fitsp.exitflag = exitflag;
fitsp.output   = output;

end