function fitsp = fitspec(obj, varargin)
% fits spin wave spectra to experimental spectral data
%
% fitsp = FITSPEC(obj, 'Option1', Value1, ...)
%
% The function uses a heuristic method to fit spin wave spectrum using a
% few simple rules to define the R-value of the fit:
%   1 All calculated spin wave modes that are outside of the measured
%     energy range will be omitted.
%   2 Spin wave modes that are closer to each other than the given energy
%     bin will be binned together and considered as one mode in the fit.
%   3 If the number of calculated spin wave modes after applying rule 1&2 
%     is larger than the observed number, the weakes simulated modes will
%     be removed from the fit.
%   4 If the number of observed spin wave modes is larger than the observed
%     number, fake spin wave modes are added with energy equal to the
%     limits of the scan; at the upper or lower limit depending on which is
%     closer to the observed spin wave mode.
% After these rules the number of observed and simulated spin wave modes
% will be equal. The R-value is defined as:
%
%       R = sqrt( nE^(-1) * sum_i_q (E_i_q_sim - E_i_q_meas)^2/sigma_i_q^2 ),
%
% where (i,q) indexing the spin wave mode and momentum. E_sim and E_meas
% are the simulated and measured spin wave energies, sigma is the standard
% deviation of the measured spin wave energy determined previously by
% fitting the inelastic peak. nE is the number of energies to fit.
%
% The R value is optimized using particle swarm algorithm to find the
% global minimum.
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
% optimizer String that determines the optimizer to use, possible values:
%               'pso'       Particle-swarm optimizer, see ndbase.pso,
%                           default.
%               'simplex'   Matlab built-in simplex optimizer, see
%                           fminsearch.
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
% imagChk   Checks that the imaginary part of the spin wave dispersion is
%           smaller than the energy bin size. Default is true.
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
% TolX          Minimum change of x when convergence reached, default
%               value is 1e-4.
% TolFun        Minimum change of the R value when convergence reached,
%               default value is 1e-5.
% MaxFunEvals   Maximum number of function evaluations, default value is
%               1e7.
% MaxIter       Maximum number of iterations for the ndbse.pso optimizer.
%               Default value is 20.
%
% Output:
%
% Output fitsp is struct type with the following fields:
% obj       Copy of the input spinw class object, with the best fitted
%           Hamiltonian.
% x         Final values of the fitted parameters, dimensions are
%           [nRun nPar]. The rows of x are sorted according to increasing R
%           values.
% R         R-value, goodness of the fit, dimensions are [nRun 1], sorted
%           in increasing order.
% exitflag  Exit flag of the fminsearch command.
% output    Output of the fminsearch command.
%
% Any option used by SPINW.SPINWAVE function are also accepted.
%
% See also SPINW.SPINWAVE, SPINW.MATPARSER, SW_EGRID, SW_NEUTRON, SW_READSPEC.
%

tid0 = swpref.getpref('tid',[]);

inpForm.fname  = {'epsilon' 'datapath' 'xmin'   'xmax'  'x0'    'func' 'plot'};
inpForm.defval = {1e-5      ' '        []       []      []      []     true  };
inpForm.size   = {[1 1]     [1 -1]     [1 -3]   [1 -4]  [1 -5]  [1 1]  [1 1] };
inpForm.soft   = {true      false      true     true    true    false  false };

inpForm.fname  = [inpForm.fname  {'tolx' 'tolfun' 'maxfunevals' 'Evect' 'nRun'}];
inpForm.defval = [inpForm.defval {1e-4   1e-5     1e3           []      1     }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]    [1 1]         [1 -6]  [1 1] }];
inpForm.soft   = [inpForm.soft   {false  false    false         false   false }];

inpForm.fname  = [inpForm.fname  {'nMax' 'hermit' 'iFact' 'lShift' 'optimizer'}];
inpForm.defval = [inpForm.defval {1e3    true     1/30     0       'pso'      }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]    [1 1]    [1 1]   [1 -7]     }];
inpForm.soft   = [inpForm.soft   {false  false    false    false   false      }];

inpForm.fname  = [inpForm.fname  {'maxiter' 'sw'  'optmem' 'tid' 'imagChk'}];
inpForm.defval = [inpForm.defval {20        1     0        tid0   true    }];
inpForm.size   = [inpForm.size   {[1 1]  	[1 1] [1 1]    [1 1]  [1 1]   }];
inpForm.soft   = [inpForm.soft   {false  	false false    false  false   }];

param = sw_readparam(inpForm, varargin{:});

% number of parameters (length of x)
nPar = max(max(length(param.xmin),length(param.xmax)),length(param.x0));
nRun = param.nRun;

% elseif isempty(param.xmax) && isempty(param.xmin)
%     x0 = rand(nRun,nPar);
% elseif isempty(param.xmax)
%     x0 = repmat(param.xmin,[nRun 1]) + rand(nRun,nPar);
% elseif isempty(param.xmin)
%     x0 = repmat(param.xmax,[nRun 1]) - rand(nRun,nPar);
% else
%     x0 = bsxfun(@times,rand(nRun,nPar),(param.xmax-param.xmin))+repmat(param.xmin,[nRun 1]);
% end

% Read experimental data
data = sw_readspec(param.datapath);

% if any sigma is zero make it 1
for ii = 1:numel(data)
    data{ii}.sigma(data{ii}.sigma==0) = 1;
end

% setup parameters for spinw.spinwave()
param0      = param;
param0.plot = false;
param0.tid  = 0;

x = zeros(nRun,nPar);
R = zeros(nRun,1);


% convert data into standard xye format
dat = struct('x',[],'y',[],'e',[]);
for ii = 1:numel(data)
    sel   = data{ii}.I~=0 & data{ii}.E~=0;
    dat.y = [dat.y;data{ii}.E(sel)];
    dat.e = [dat.e;data{ii}.sigma(sel)];
end
dat.x = (1:numel(dat.y))';
dat.y = dat.y(:);
dat.e = dat.e(:);

sw_status(0,1,param.tid,'Fitting spin wave spectra');

idx = 1;
idxAll = 1;

fidSave = obj.fileid;
obj.fileid (0);

while idx <= nRun
    %try
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
    
    switch param.optimizer
        case 'pso'
            [x(idx,:),fVal, output(idx)] = ndbase.pso(dat,@(x,p)spec_fitfun(obj, data, param.func, p, param0),x0,'lb',param.xmin,'ub',param.xmax,...
                'TolX',param.tolx,'TolFun',param.tolfun,'MaxFunEvals',param.maxfunevals,'MaxIter',param.maxiter);
            R(idx) = sqrt(sum(((dat.y-fVal(:))./dat.e).^2)/numel(fVal));
            
        case 'simplex'
            [x(idx,:), R(idx), ~, output(idx)] = sw_fminsearchbnd(@(p)spec_fitfun(obj, data, param.func, p, param0),x0,param.xmin,param.xmax,...
                optimset('TolX',param.tolx,'TolFun',param.tolfun,'MaxFunEvals',param.maxfunevals,'Display','off'));
        otherwise
            error('spinw:fitspec:WrongOption','The given optimizer is not supported!')
            
            %case 'both'
            %    [xPSO, fPSO] = ndbase.pso(dat,@(x,p)spec_fitfun(obj, data, param.func, p, param0),x0,'lb',param.xmin,'ub',param.xmax,...
            %        'TolX',param.tolx,'TolFun',param.tolfun,'MaxFunEvals',param.maxfunevals,'MaxIter',20);
            %    [x(idx,:),fVal, output(ii)] = ndbase.lm3(dat,@(x,p)spec_fitfun(obj, data, param.func, p, param0),xPSO,'lb',param.xmin,'ub',param.xmax,...
            %        'TolX',param.tolx,'TolFun',param.tolfun,'MaxFunEvals',param.maxfunevals);
            
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
    sw_status(idx/nRun*100,0,param.tid);
end

sw_status(100,2,param.tid);

% Sort results
[R, sortIdx] = sort(R);
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

obj.fileid(fidSave);

% Store all output in a struct variable.
fitsp.obj      = copy(obj);
fitsp.x        = x;
fitsp.R        = R;
%fitsp.exitflag = exitflag;
fitsp.output   = output;

end

function [yCalc, pHandleOut] = spec_fitfun(obj, dataCell, parfunc, x, param)
% calculates the spin wave energies to directly compare to data
%
% [yCalc, pHandle] = SW_FITFUN(obj, data, param, swfunc, x)
%
% swobj         spinw type object contains the magnetic structure and the
%               Hamiltonian for the spin wave calculation.
% data          Structure contains the experimental data, read by
%               sw_readspec.
% param         Paramters for the spinwave function.
% swfunc        Function to change the Hamiltonian in obj, has the
%               following form:
%                   obj = swfunc(obj,x);
% x             Actual parameter values.
%

param.nPoints = 50;


parfunc(obj,x);

% Number of different correlation functions measured.
nConv = numel(dataCell);
yCalc = [];
R     = 0;

sim = struct;
Qc = 0;

for ii = 1:nConv
    % Select one type of correlation
    data = dataCell{ii};
    
    % calculate spin-spin correlation function
    spec = obj.spinwave(data.Q,'fitmode',true,'hermit',param.hermit,...
        'tid',0,'optMem',param.optmem,'tid',param.tid,'gtensor',any(obj.single_ion.g));
    % calculate neutron scattering cross section
    spec = sw_neutron(spec,'n',data.n,'pol',data.corr.type{1}(1) > 1);
    % bin the data along energy
    spec = sw_egrid(spec,'component',data.corr,'Evect',param.Evect,'imagChk',param.imagChk);
    % generate center bin
    spec.Evect = (spec.Evect(1:(end-1))+spec.Evect(2:end))/2;
    
    nQ = size(data.Q,2);
    pHandle = [];
    
    for jj = 1:nQ
        Erange    = (spec.Evect >= data.minE(jj)) & (spec.Evect <= data.maxE(jj));
        Evect     = spec.Evect(Erange);
        swConv    = spec.swConv(Erange,jj);
        
        data.nii  = data.nMode(jj);
        
        [~, idx]  = sort(swConv,'descend');
        sim.nMode = min(data.nii,sum(swConv>0));
        idx       = idx(1:sim.nMode);
        sim.I     = swConv(idx);
        sim.E     = Evect(idx);
        % sort calculated magnon energies in increasing order
        [sim.E, idx2] = sort(sim.E);
        sim.I     = sim.I(idx2);
        
        data.Eii  = data.E(1:data.nii,jj)';
        data.Iii  = data.I(1:data.nii,jj)';
        data.sii  = data.sigma(1:data.nii,jj)';
        data.wii  = data.sii.^(-2);
        
        simVect   = [repmat(data.minE(jj),[1 data.nii]) sim.E repmat(data.maxE(jj),[1 data.nii])];
        
        nSlip     = (data.nii+sim.nMode+1);
        d = zeros(1,nSlip);
        for kk = 1:nSlip
            % R-value calculated as weight(Data)*(E(simulation)-E(Data))^2
            d(kk) = sum(data.wii.*(simVect((1:data.nii)+kk-1)-data.Eii).^2);
        end
        [mind,Eidx] = min(d);
        yCalc = [yCalc;simVect((1:data.nii)+Eidx-1)'];
        R = R + mind;
        
        if param.plot
            plot(jj+Qc+sim.E*0,sim.E,'ko');
            hold on
            for kk = 1:length(sim.E)
                cPoints = sw_circle([jj+Qc sim.E(kk) 0]',[0 0 1]',sqrt(sim.I(kk))*param.iFact,param.nPoints);
                pHandle(end+1) = plot(cPoints(1,:),cPoints(2,:)); %#ok<*AGROW>
                set(pHandle(end),'Color','k');
            end
            
            pHandle(end+1) = plot(jj+Qc+data.Eii*0,data.Eii,'r.');
            for kk = 1:length(data.Eii)
                cPoints = sw_circle([jj+Qc data.Eii(kk) 0]',[0 0 1]',sqrt(data.Eii(kk))*param.iFact,param.nPoints);
                pHandle(end+1) = plot(cPoints(1,:),cPoints(2,:));
                set(pHandle(end),'Color','r');
                pHandle(end+1) = line([0 0]+jj+Qc,[-1 1]*data.sii(kk)*param.sw+data.Eii(kk),'Color','r','LineWidth',2);
            end
            
            pHandle(end+1) = line([-0.2 0.2]+jj+Qc,data.minE(jj)*[1 1],'Color','g');
            pHandle(end+1) = line([-0.2 0.2]+jj+Qc,data.maxE(jj)*[1 1],'Color','g');
            qLabel = sprintf('(%4.2f,%5.2f,%5.2f)',spec.hkl(:,jj));
            pHandle(end+1) = text(jj+Qc-0.1,(param.Evect(end)-param.Evect(1))*0.7+param.Evect(1)+param.lShift,qLabel,'rotation',90,'fontsize',12);
        end
    end
    Qc = Qc + nQ;
end

R = R/numel(yCalc);

if param.plot
    text(0.05,0.9,['x = [' sprintf('%6.4f ',x) sprintf(']\nRw = %6.4f',sqrt(R))],'Units','normalized','fontsize',12);
    axis([0.5 Qc+0.5 param.Evect(1) param.Evect(end)]);
    legend(pHandle(1:2),'simulation','data')
    xlabel('Scan index')
    ylabel('Energy transfer (meV)')
    title('Spin wave dispersion fit')
    set(gca,'XTick',1:nQ)
    hold off
    drawnow;
end

if nargout>1
    pHandleOut = pHandle;
end

end