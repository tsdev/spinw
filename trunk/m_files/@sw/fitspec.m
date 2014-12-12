function fitsp = fitspec(obj, varargin)
% fits spin wave spectra to experimental spectral data
%
% fitsp = FITSPEC(obj, 'Option1', Value1, ...)
%
% Options:
%
% epsilon   Small number that controls wether the magnetic structure is
%           incommensurate or commensurate, default value is 1e-5.
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
% plot      If true, the measured dispersion is plotted together with the
%           fit. Default is true.
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
% See also SW.SPINWAVE, SW_EGRID, SW_NEUTRON, SW_READSPEC, FMINSEARCH.
%

inpForm.fname  = {'epsilon' 'datapath' 'xmin'   'xmax'  'x0'    'func' 'plot'};
inpForm.defval = {1e-5      ' '        []       []      []      []     true  };
inpForm.size   = {[1 1]     [1 -1]     [1 -3]   [1 -4]  [1 -5]  [1 1]  [1 1] };
inpForm.soft   = {1         0          1        1       1       0      0     };

inpForm.fname  = [inpForm.fname  {'tolx' 'tolfun' 'maxfunevals' 'Evect' 'nRun'}];
inpForm.defval = [inpForm.defval {1e-4   1e-5     1e7           []      1     }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]    [1 1]         [1 -6]  [1 1] }];
inpForm.soft   = [inpForm.soft   {0      0        0             0       0     }];

inpForm.fname  = [inpForm.fname  {'nMax' 'hermit'}];
inpForm.defval = [inpForm.defval {1e3    true    }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]   }];
inpForm.soft   = [inpForm.soft   {0      0       }];

param = sw_readparam(inpForm, varargin{:});

% number of parameters (length of x)
nPar = max(max(length(param.xmin),length(param.xmax)),length(param.x0));
nRun = param.nRun;

% Initial parameters are random if param.x0 is undefined.
if ~isempty(param.x0)
    nRun = 1;
end
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

param0      = param;
param0.plot = false;

x = zeros(nRun,nPar);
R = zeros(nRun,1);
%exitflag = zeros(nRun,1);
%output   = struct;

sw_status(0,1);

idx = 1;
idxAll = 1;

fidSave = obj.fileid;
obj.fileid (0);

while idx <= nRun
    try
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
        
        [x(idx,:), R(idx), exitflag(idx), output(idx)] = sw_fminsearchbnd(@(x)sw_fitfun(obj, data, param.func, x, param0),x0,param.xmin,param.xmax,...
            optimset('TolX',param.tolx,'TolFun',param.tolfun,'MaxFunEvals',param.maxfunevals,'Display','off'));

        idx = idx + 1;
    catch
        warning('Hamiltonian is not compatible with the magnetic ground state!');
        if ~isempty(param.x0)
            error('sw:fitspec:WrongInput','x0 input incompatible with ground state!');
        end
        
        if idxAll >= param.nMax
            warning('sw:fitspec:WrongInput','Maximum number of wrong runs (param.nMax) reached, calculation is terminated!');
            break;
        end
    end
    idxAll = idxAll + 1;
    sw_status(idx/nRun*100);
end
sw_status(100,2);

% Sort results
[R, sortIdx] = sort(R);
x = x(sortIdx,:);

% Draw plot of the final result if requested
if param.plot
    figure;
    param0.plot = true;
    % plot the best result if requested
    sw_fitfun(obj, data, param.func, x(1,:), param0);
end

% set the best fit to the sw object
obj = param.func(obj,x(1,:));

obj.fileid(fidSave);

% Store all output in a struct variable.
fitsp.obj      = copy(obj);
fitsp.x        = x;
fitsp.R        = R;
fitsp.exitflag = exitflag;
fitsp.output   = output;

end

function [R, pHandle] = sw_fitfun(obj, dataCell, parfunc, x, param)
% [R, pHandle] = SW_FITFUN(obj, data, param, swfunc, x) calculates the
% agreement factor between simulated and measured data.
%
% swobj         sw type object contains the magnetic structure and the
%               Hamiltonian for the spin wave calculation.
% data          Structure contains the experimental data, read by
%               sw_readspec.
% param         Paramters for the spinwave function.
% swfunc        Function to change the Hamiltonian in obj, has the
%               following form:
%                   obj = swfunc(obj,x);
% x             Actual parameter values.
%

param.iFact   = 1/30;
param.nPoints = 50;

obj = parfunc(obj,x);

%nExt = double(obj.mag_str.N_ext);

% Number of different correlation functions measured.
nConv = numel(dataCell);
R  = 0;

sim = struct;
Qc = 0;

for ii = 1:nConv
    % Select one type of correlation
    data = dataCell{ii};
    
    % calculate spin-spin correlation function
    spec = obj.spinwave(data.Q,'fitmode',true,'hermit',param.hermit);
    % calculate neutron scattering cross section
    spec = sw_neutron(spec,'n',data.n,'pol',data.corr.type{1}(1) > 1);
    % bin the data along energy
    spec = sw_egrid(spec,'component',data.corr,'Evect',param.Evect);
    
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
        data.wii  = data.w(1:data.nii,jj)';
        
        simVect   = [repmat(data.minE(jj),[1 data.nii]) sim.E repmat(data.maxE(jj),[1 data.nii])];
        
        nSlip     = (data.nii+sim.nMode+1);
        d = zeros(1,nSlip);
        for kk = 1:nSlip
            % R-value calculated as width(Data)*(E(simulation)-E(Data))^2
            d(kk) = sum(data.wii.*(simVect((1:data.nii)+kk-1)-data.Eii).^2);
        end
        R = R + min(d);
        
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
                pHandle(end+1) = line([0 0]+jj+Qc,[-data.wii(kk)/2 data.wii(kk)/2]+data.Eii(kk),'Color','r','LineWidth',2);
            end
            
            pHandle(end+1) = line([-0.2 0.2]+jj+Qc,data.minE(jj)*[1 1],'Color','g');
            pHandle(end+1) = line([-0.2 0.2]+jj+Qc,data.maxE(jj)*[1 1],'Color','g');
            qLabel = sprintf('(%4.2f,%5.2f,%5.2f)',spec.hkl(:,jj));
            pHandle(end+1) = text(jj+Qc-0.1,(param.Evect(end)-param.Evect(1))*0.7+param.Evect(1),qLabel,'rotation',90,'fontsize',12);
        end
    end
    Qc = Qc + nQ;
end

if param.plot
    text(0.05,0.9,['x = [' sprintf('%6.4f ',x) sprintf(']\nRw = %6.4f',R)],'Units','normalized','fontsize',12);
    axis([0.5 Qc+0.5 param.Evect(1) param.Evect(end)]);
    legend(pHandle(1:2),'simulation','data')
    xlabel('Scan index')
    ylabel('Energy transfer (meV)')
    title('Spin wave dispersion fit')
    set(gca,'XTick',1:nQ)
    hold off
    drawnow;
end

end