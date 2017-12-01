function varargout = powder_optimiser(obj, data, Q_slices, p0, varargin)

%% Function workflow

%%% Preprocess the data - should the user do this before?

% Remove the elastic line if it is there.
% We set all outside range to NaN.

%%% We chooose N q (NQ) slices to optimise on.
% Given in q_slices.
% `q_slices = {[Q1_start, Q1_end, nQ1], [Q2_start, Q2_end, nQ2], ...}`

%%% Simulation options.
% Each Q slice will be of P points (nQ1, nQ2, .... ).

% We start the pso optimiser with the current values. Evaluate for
% q_slices.

%%%  PSO Loop
% In each loop we will use the fibonaci option and have nRand points.

% The computed sprectum will be convolved and compared to the data. R^2?

% At the update stage `matparser` will be used to update the exchanges.

%%% LSQ Loop
% When a solution is approximated a least squares method will be employed.

%%% Return.
% The optimised object will be returned with an R value.


%% Input parameters.

% Initial spinW object

% Data to be optimised (single dat.x, dat.y, dat.z, dat.e for now)

% Q Points - `q_slices = {[Q1_start, Q1_end, nQ1], [Q2_start, Q2_end, nQ2], ...};`

% OPTIONAL

% Exchange labels to be fitted. `exch_lab = {'J1', 'J2', 'J3(1,2)' ...};`

% Starting points. Do not use the current obj.

% PSO options.

%% Start of Code

pref = swpref;

title0 = 'Powder LSWT spectrum';
tid0   = pref.tid;

% These are the powspec options.
inpForm.fname  = {'nRand' 'T'   'formfact' 'formfactfun' 'tid' 'nInt'  'dE'    'exch_lab'};
inpForm.defval = {100     T0    false      @sw_mff       tid0  1e3     2       obj.matrix.label};
inpForm.size   = {[1 1]   [1 1] [1 -2]     [1 1]         [1 1] [1 1]   [1, 1]  [1, -6]};

inpForm.fname  = [inpForm.fname  {'hermit' 'gtensor' 'title' 'specfun' 'imagChk'}];
inpForm.defval = [inpForm.defval {true     false     title0  @spinwave  true    }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]     [1 -3]  [1 1]      [1 1]   }];

inpForm.fname  = [inpForm.fname  {'extrap' 'fibo' 'optmem' 'binType' 'component'}];
inpForm.defval = [inpForm.defval {false    false  0        'ebin'    'Sperp'    }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]  [1 1]    [1 -4]     [1 -5]    }];

inpForm.fname  = [inpForm.fname  {'fid'}];
inpForm.defval = [inpForm.defval {-1   }];
inpForm.size   = [inpForm.size   {[1 1]}];

warning('off','sw_readparam:UnreadInput')
param_PS  = sw_readparam(inpForm, varargin{:});
warning('on','sw_readparam:UnreadInput')


%TODO: define these vars properly...
exch_lab = param_PS.exch_lab;
Np = length(exch_lab);
oNp = ones(1,Np);

% PSO options.
inpForm.fname  = {'Display' 'TolFun' 'TolX' 'MaxIter' 'SwarmC1' 'seed'        };
inpForm.defval = {'off'     1e-3     1e-3   100*Np    2.8      sum(100*clock) };
inpForm.size   = {[1 -1]    [1 1]    [1 1]  [1 1]     [1 1]    [1 1]          };

inpForm.fname  = [inpForm.fname  {'SwarmC2' 'PopulationSize' 'lb'      'ub'     'autoTune' 'k0'}];
inpForm.defval = [inpForm.defval {1.3       25               -1e5*oNp   1e5*oNp false      1   }];
inpForm.size   = [inpForm.size   {[1 1]     [1 1]            [1 Np]     [1 Np]  [1 1]     [1 1]}];

warning('off','sw_readparam:UnreadInput')
param_PSO  = sw_readparam(inpForm, varargin{:});
warning('on','sw_readparam:UnreadInput')


%% Slice the experimental data

data_Q = data.x;
data_E = data.y;
data_I = data.z;
data_V = date.e;

% This is used for the simulaiton.
data_S = struct('Q',[],'E',[],'I',[],'V',[],'ind',[],'Q_edge',[],'Q_edge_ind',[]);

for i = 1:length(Q_slices)
    
    % Split up the data.
    Q_slice = Q_slices{i};
    this_Q = data_Q((data_Q >= Q_slice(1)) & (data_Q <= Q_slice(2)));
    this_E = data_E((data_Q >= Q_slice(1)) & (data_Q <= Q_slice(2)));
    this_I = data_I((data_Q >= Q_slice(1)) & (data_Q <= Q_slice(2)));
    this_V = data_V((data_Q >= Q_slice(1)) & (data_Q <= Q_slice(2)));
    
    % Do interpolation.
    u_q = this_Q;
    n_q = length(u_q);
    if n_q == 0
        warning('spinw:powder_optimiser:InvalidQRange','The Q range has no q points.')
        continue
    elseif n_q > Q_slice(3)
        warning('spinw:powder_optimiser:InvalidQRange','The Q range has too many q points.')
    else
        Q_slice(3) = n_Q;
    end
    
    [N, Q_EDGES, BIN] = histcounts(this_Q, Q_slice(3)); %#ok<ASGLU>
    this_E = accumarray(BIN(:),this_E(:),[],@mean);
    this_I = accumarray(BIN(:),this_I(:),[],@mean);
    
    N = accumarray(BIN(:),ones(size(this_Q)),[],@sum,NaN);
    this_V = accumarray(BIN(:),this_V(:),[],@norm)./N(:);
    this_Q = accumarray(BIN(:),this_Q(:),[],@mean);
    
    data_S.Q = [data_S.Q; this_Q(:)];
    data_S.Q_edge = [data_S.Q_edge; Q_EDGES(:)];
    data_S.Q_edge_ind = [data_S.Q_edge_ind; i*ones(length(Q_EDGES),1)];
    data_S.E = [data_S.E; this_E(:)];
    data_S.I = [data_S.I; this_I(:)];
    data_S.V = [data_S.V; this_V(:)];
    data_S.ind = [data_S.ind; i*ones(length(this_Q),1)];
    Q_slices{i} = Q_slice;
end

dat = struct('x',zeros(size(data_S)),'y',data_S.I,'e',data_S.V);

p0 = [p0(:); param_PS.dE];
[pOpt, fVal, stat] = ndbase.pso(dat,@do_simulation_loop,p0,param_PSO);

obj_final = obj.copy;
obj_final.matparser(exch_lab,pOpt);

varargout{1} = obj_final;
if nargout > 1
    this_param = rmfield(param_PS,{'dE', 'exch_lab'});
    final_spectra = obj_final.powspec(data_Q, 'Evect', data_E, this_param);
    final_spectra = sw_instrument(final_spectra,'dE', param_PS.dE);
    varargout{2} =  sum(((data_I(:) - final_spectra.swConv(:))./data_V(:)).^2)/(length(data_I(:)) - Np);
end
if nargout == 3
    varargout{2} = stat;
end


    function y_data = do_simulation_loop(~, p)
        
        dE = p(end);
        p = p(1:end-1);
        
        obj_local = obj.copy;
        obj_local.matparser(exch_lab,p);
        
        runs = unique(data_S.ind);
        y_data = [];
        for j = 1:runs
            param_local = rmfield(param_PS,{'dE', 'exch_lab'});
            this_spectra = obj_local.powspec(data_S.Q_edge(data_S.Q_edge_ind == j),...
                'Evect',data_S.Q_edge(data_S.Q_edge_ind == j),param_local);
            % We need to do an instrumental convolution at this points....
            this_spectra = sw_instrument(this_spectra,'dE', dE);
            y_data = [y_data(:); this_spectra.swConv];
        end
    end
end