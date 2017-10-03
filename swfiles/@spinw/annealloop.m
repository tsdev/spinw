function stat = annealloop(obj, varargin)
% parameter sweep for simulated annealing
% 
% ### Syntax
% 
% `stat = anneal(obj,Name,Value)`
% 
% ### Description
% 
% `stat = annealloop(obj,Name,Value)` performs simulated annealing while
% stepwise changing a selected parameter such as temperature or magnetic
% field while measuring thermodynamic properties. The function has the same
% parameters as [spinw.anneal] plus an additional
%  
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Name-Value Pair Arguments
% 
% `'func'`
% : Function that changes the parameters in the spinw object in every
%   loop. Default function is to change the temperature:
%   ```
%   @temperature(obj,x)
%   ```
%   The function takes two inputs a [spinw] object and a parameter value
%   (ir values in a vector) and changes the correspondign property inside
%   the object.
% 
% `'x'`
% : Matrix of values of the loop parameter, with dimensions of
%   $[n_{par}\times n_{loop}]$. Default value is 1. In the i-th loop the
%   loop function is called as:
%   ```
%   func(obj,x(:,i));
%   ```
% 
% `'saveObj'`
% : If `true`, the spinw object is saved after every annealing step for
%   debugging purposes. Default is `false`.
%
% `'tid'`
% : Determines if the elapsed and required time for the calculation is
%   displayed. The default value is determined by the `tid` preference
%   stored in [swpref]. The following values are allowed (for more details
%   seee [sw_status]):
%   * `0` No timing is executed.
%   * `1` Display the timing in the Command Window.
%   * `2` Show the timing in a separat pup-up window.
% 
% ### Output Arguments
% 
% Same output as of [spinw.anneal], just the struct is packaged into a cell
% with $n_{loop}$ number of elements.
%
% ### Reference
%
%    Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by
%    Simulated Annealing. _Science, 220_, 671-680.
% 
% ### See Also
% 
% [spinw] \| [spinw.anneal] \| [sw_fsub] \| [sw_fstat]

func0 = @(obj,T)obj.temperature(T);

nExt   = double(obj.magstr.N_ext);

inpForm.fname  = {'initT' 'endT' 'cool'        'nMC' 'nStat' };
inpForm.defval = {100     1e-2   @(T) (0.92*T) 100   100     };
inpForm.size   = {[1 1]   [1 1]  [1 1]         [1 1] [1 1]   };
inpForm.soft   = {false   false  false         false  false  };

inpForm.fname  = [inpForm.fname  {'spinDim' 'nORel' 'nExt' 'subLat'     }];
inpForm.defval = [inpForm.defval {3         0       nExt   []           }];
inpForm.size   = [inpForm.size   {[1 1]     [1 1]   [1 3]  [1 -1]       }];
inpForm.soft   = [inpForm.soft   {false     false   false  true         }];

inpForm.fname  = [inpForm.fname  {'fStat'   'fSub'   'random' 'title'   }];
inpForm.defval = [inpForm.defval {@sw_fstat @sw_fsub false    ''        }];
inpForm.size   = [inpForm.size   {[1 1]     [1 1]    [1 1]    [1 -2]	}];
inpForm.soft   = [inpForm.soft   {false     false    false    true      }];

inpForm.fname  = [inpForm.fname  {'fineT' 'rate' 'boundary'         }];
inpForm.defval = [inpForm.defval {0       0.1    {'per' 'per' 'per'}}];
inpForm.size   = [inpForm.size   {[1 1]   [1 1]  [1 3]              }];
inpForm.soft   = [inpForm.soft   {false   false  false              }];

inpForm.fname  = [inpForm.fname  {'x'     'func' 'saveObj' 'tid'}];
inpForm.defval = [inpForm.defval {1       func0  false     -1   }];
inpForm.size   = [inpForm.size   {[-3 -4] [1 1]  [1 1]     [1 1]}];
inpForm.soft   = [inpForm.soft   {false   false  false     false}];

param = sw_readparam(inpForm,varargin{:});

if param.tid == -1
    param.tid = swpref.getpref('tid',[]);
end

% check output
param.verbosity = 0;

nLoop = size(param.x,2);

aRes = cell(1,nLoop);

sw_status(0,1,param.tid,'Parameter sweep for simulated annealing');

stat = struct;

for ii = 1:nLoop
    
    % change parameters
    param.func(obj,param.x(:,ii));
    param.initT = obj.temperature;
    param.endT  = obj.temperature;
    
    % do the annealing procedure
    warnState = warning('off','sw_readparam:UnreadInput');
    paramF = param;
    paramF.fastmode = ~param.saveObj;
    aRes{ii} = obj.anneal(paramF);
    warning(warnState);
    
    % save sublattice list
    param.subLat = aRes{ii}.param.subLat;
    
    if strcmp(func2str(param.fStat),'sw_fstat')
        if ii>1
            % avgM dimensions are [3 nMagExt]
            stat.avgM = cat(3,stat.avgM,aRes{ii}.avgM);
            % stdM dimensions are [3 nMagExt]
            stat.stdM = cat(3,stat.stdM,aRes{ii}.stdM);
            % avgE average system energy per spin over nStat runs, scalar.
            stat.avgE = [stat.avgE aRes{ii}.avgE];
            % stdE standard deviation of the system energy per spin
            stat.stdE = [stat.stdE aRes{ii}.stdE];
            % Cp heat capacity of the sample: (<E^2>-<E>^2)/kB/T^2.
            stat.Cp = [stat.Cp aRes{ii}.Cp];
            % Chi magnetic susceptibility of the sample: (<M^2>-<M>^2)/kB/T.
            stat.Chi = [stat.Chi aRes{ii}.Chi];
        else
            stat.avgM = aRes{ii}.avgM;
            stat.stdM = aRes{ii}.stdM;
            stat.avgE = aRes{ii}.avgE;
            stat.stdE = aRes{ii}.stdE;
            stat.Cp   = aRes{ii}.Cp;
            stat.Chi  = aRes{ii}.Chi;
        end
        
    else
        % TODO
        warning('TODO');
    end
    
    sw_status(ii/nLoop*100,0,param.tid);
end

% save extra information
stat.obj = copy(obj);
stat.state = aRes;
stat.param = param;

sw_status(100,2,param.tid);

end
