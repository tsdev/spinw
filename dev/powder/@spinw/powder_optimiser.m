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

%%
% NOTE p0 is [vals for matparser; IntensityFactor];
%

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
T0 = 0;

% These are the powspec options.
inpForm.fname  = {'nRand' 'T'   'formfact' 'formfactfun' 'tid' 'nInt'  'dE'    'exch_lab'};
inpForm.defval = {100     T0    false      @sw_mff       0  1e3     2       obj.matrix.label};
inpForm.size   = {[1 1]   [1 1] [1 -2]     [1 1]         [1 1] [1 1]   [1, 1]  [1, -6]};

inpForm.fname  = [inpForm.fname  {'hermit' 'gtensor' 'title' 'specfun' 'imagChk'}];
inpForm.defval = [inpForm.defval {true     false     title0  @spinwave  true    }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]     [1 -3]  [1 1]      [1 1]   }];

inpForm.fname  = [inpForm.fname  {'extrap' 'fibo' 'optmem' 'binType' 'component'}];
inpForm.defval = [inpForm.defval {false    true   0        'cbin'    'Sperp'    }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]  [1 1]    [1 -4]     [1 -5]    }];

inpForm.fname  = [inpForm.fname  {'fid'}];
inpForm.defval = [inpForm.defval {0   }];
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
if isempty(data.e)
    data_V = 1./data.z;
else
    data_V = date.e;
end
% This is used for the simulaiton.
data_S = struct('Q',[],'E',[],'I',[],'V',[]);

k = 1;
for i = 1:length(Q_slices)
    
    % Split up the data.
    Q_slice = Q_slices{i};
    q_ind = (data_Q(1,:) >= Q_slice(1)) & (data_Q(1,:) <= Q_slice(2));
    this_Q = data_Q(:,q_ind);
    
    this_E = data_E(:,q_ind);
    this_I = data_I(:,q_ind);
    this_V = data_V(:,q_ind);
    
    % Do interpolation.
    n_q = sum(q_ind);
    if n_q == 0
        warning('spinw:powder_optimiser:InvalidQRange','The Q range has no q points.')
        % We skip this itteration. i*ones doesn't matter as we une unique
        % indexing later....
        continue
    elseif n_q > Q_slice(3)
        % We will be doing binning.
        warning('spinw:powder_optimiser:InvalidQRange','The Q range has too many q points. Binning to %f points.',Q_slice(3))
    else
        % We will be doing binning, but we shouldnt be doing binning....
        Q_slice(3) = n_q;
    end
    
    [N, Q_EDGES, BIN] = histcounts(this_Q, Q_slice(3)); %#ok<ASGLU>
    this_E = accumarray([BIN(:) reshape(repmat((1:size(BIN,1))',1,size(BIN,2)),[],1)],this_E(:),[],@mean)';
    this_I = accumarray([BIN(:) reshape(repmat((1:size(BIN,1))',1,size(BIN,2)),[],1)],this_I(:),[],@mean)';
    
    N = accumarray(BIN(:),reshape(ones(size(this_Q)),1,[]),[],@sum,NaN);
    this_V = bsxfun(@rdivide,accumarray([BIN(:) reshape(repmat((1:size(BIN,1))',1,size(BIN,2)),[],1)],this_V(:),[],@norm)',N');
%     this_Q = accumarray([BIN(:) reshape(repmat((1:size(BIN,1))',1,size(BIN,2)),[],1)],this_Q(:),[],@mean)';
    
    data_S(k).Q = Q_EDGES(1:end-1) + diff(Q_EDGES);
    data_S(k).E = this_E;
    data_S(k).I = this_I;
    data_S(k).V = this_V;
    Q_slices{i} = Q_slice;
    k = k+1;
end

dat = struct('x',zeros(numel([data_S.I]),1),'y',reshape([data_S.I],[],1),'e',reshape([data_S.V],[],1));
dat.x = dat.x(~isinf(dat.e));
dat.y = dat.y(~isinf(dat.e));
dat.e = dat.e(~isinf(dat.e));

p0 = [p0(:); param_PS.dE];
param_PSO.lb = [param_PSO.lb 0 0];
param_PSO.ub = [param_PSO.ub p0(3:4)'*100];

[pOpt, fVal, stat] = ndbase.pso(dat,@do_simulation_loop,p0,param_PSO); %#ok<ASGLU>
[pOpt, fVal, stat] = ndbase.lm3(dat,@do_simulation_loop,pOpt,'ub',param_PSO.ub,'lb',param_PSO.lb); %#ok<ASGLU>

obj_final = obj.copy;
obj_final.matparser('param',pOpt(1:end-2),'mat',exch_lab);

varargout{1} = obj_final;
varargout{2} = pOpt;

% if nargout > 1
%     this_param = rmfield(param_PS,{'dE', 'exch_lab'});
%     this_param.Evect = data.E(:);
%     final_spectra = obj_final.powspec(data_Q(:), this_param);
%     final_spectra = sw_instrument(final_spectra,'dE', pOpt(end));
%     final_spectra.swConv = final_spectra.swConv * pOpt(end-1);
% end
if nargout == 3
    varargout{3} = stat;
end


    function y_data = do_simulation_loop(dummy, p) %#ok<INUSL>
        
        dE = p(end);
        I  = p(end-1);
        p  = p(1:end-2);
        
        obj_local = obj.copy;
        obj_local.matparser('param',p,'mat',exch_lab);
        
        runs = length(data_S);
        y_data = [];
        param_local = rmfield(param_PS, {'dE', 'exch_lab'});

        for j = 1:runs
            for l = 1:(length(data_S(j).Q))
                % Set the Evector
                param_local.Evect = data_S(j).E(:,l)';
                % Calculate the spectra
                this_spectra = obj_local.powspec(data_S(j).Q(l), param_local);
                % Convolve
                this_spectra = sw_instrument(this_spectra, 'dE', dE,'fid',param_local.fid);
                % Add the data to the result.
                y_data = [y_data(:); this_spectra.swConv(~isinf(data_S(j).V(:,l)))];
            end
        end
        y_data = y_data*I;
    end
end