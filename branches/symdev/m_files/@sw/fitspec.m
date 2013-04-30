function fitsp = fitspec(obj, varargin)
% fits spin wave spectra to experimental data
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
% Output is struct type with the following fields:
% obj       Copy of the input sw class object, with the fitted Hamiltonian.
% x         Final values of the fitted parameters.
% R         R-value, goodness of the fit.
% exitflag  Exit flag of the fminsearch command.
% output    Output of the fminsearch command.
%
% All other options be used that SPINWAVE or SWINC functions accept. For
% the calculation of the spectra the SPINWAVE function is called, unless
% the magnetic ordering wave vector is incommensurate.
%
% See also SW.SWINC, SW.SPINWAVE, SW_CONV, SW_NEUTRON, SW_READSPEC, FMINSEARCH.
%

inpForm.fname  = {'epsilon' 'datapath' 'xmin'   'xmax'  'x0'    'func' 'plot'};
inpForm.defval = {1e-5      ' '        []       []      []      []     true  };
inpForm.size   = {[1 1]     [1 -1]     [1 -3]   [1 -4]  [1 -5]  [1 1]  [1 1] };
inpForm.soft   = {1         0          1        1       1       0      0     };

inpForm.fname  = [inpForm.fname  {'tolx' 'tolfun' 'maxfunevals' 'Evect'}];
inpForm.defval = [inpForm.defval {1e-4   1e-5     1e7           []     }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]    [1 1]         [1 -6] }];
inpForm.soft   = [inpForm.soft   {0      0        0             0      }];

param = sw_readparam(inpForm, varargin{:});

% number of parameters (length of x)
nPar = max(max(length(param.xmin),length(param.xmax)),length(param.x0));

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

% Read experimental data
data = sw_readspec(param.datapath);

param0      = param;
param0.plot = false;


[x, ~, exitflag, output] = sw_fminsearchbnd(@(x)sw_fitfun(obj, data, param.func, x, param0),x0,param.xmin,param.xmax,...
    optimset('TolX',param.tolx,'TolFun',param.tolfun,'MaxFunEvals',param.maxfunevals,'Display','off'));

% Draw plot of the final result if requested
if param.plot
    figure;
    param0.plot = true;
end

R = sw_fitfun(obj, data, param.func, x, param0);

obj = param.func(obj,x);

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
% param         Paramters for the spin wave calculation (sw_swinc or
%               sw_spinwave depending on whether the structure commensurate
%               or incommensurate).
% swfunc        Function to change the Hamiltonian in obj, has the
%               following form:
%                   obj = swfunc(obj,x);
% x             Actual parameter values.
%

param.iFact   = 1/30;
param.nPoints = 50;


obj = parfunc(obj,x);

nExt = double(obj.mag_str.N_ext);

% Choose spin wave calculation code based on the ordering wave vector.
kExt = obj.mag_str.k.*nExt;
if any(abs(kExt-round(kExt))>param.epsilon)
    swfunc = @obj.swinc;
else
    swfunc = @obj.spinwave;
end

% Number of different correlation functions measured.
nConv = numel(dataCell);
R  = 0;

sim = struct;
Qc = 0;

for ii = 1:nConv
    % Select one type of correlation
    data = dataCell{ii};
    
    % calculate spin-spin correlation function
    spec = swfunc(data.Q,'fitmode',true);
    % calculate neutron scattering cross section
    spec = sw_neutron(spec,'n',data.n,'pol',data.corr.type{1}(1) > 1);
    % bin the data along energy
    spec = sw_conv(spec,'convmode',data.corr,'Evect',param.Evect);
    
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
        
        data.Eii  = data.E(1:data.nii,jj)';
        data.Iii  = data.I(1:data.nii,jj)';
        data.wii  = data.w(1:data.nii,jj)';
        
        simVect   = [repmat(data.minE(jj),[1 data.nii]) sim.E repmat(data.maxE(jj),[1 data.nii])];
        
        nSlip     = (data.nii+sim.nMode+1);
        d = zeros(1,nSlip);
        for kk = 1:nSlip
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