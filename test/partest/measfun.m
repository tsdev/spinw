function measfun(fun,argin,usemex,nMemSlice,nRun,fName)
% measure function execution time
%
% `measfun(fun,argin,usemex,nMeMSlice,nRun,fName)` time the function
%
% `measfun` draws the header
%
% `measfun(fName)` shows the results stored in the file referenced by
% `fName`.
%

result = [];
header = sprintf('%-30s | %-15s | %-10s | %-7s | %-7s | %-7s | %-7s | %-7s | %-7s | %-7s | %-7s\n',...
    'Function name','Speed (px/sec)','Time (s)','nQ','nMode','nWorker','nThread','nSlice','nRun','useMex','hermit');

if nargin < 3
    % print header
    fprintf(header);
    % read the data
    % just print the data
    if nargin>0
        if ischar(fun)
            temp = load(fun);
            fn = fieldnames(temp);
            result = temp.(fn{1});
        else
            result = fun;
        end
    end
    if nargin>1
        % sort column
        col      = argin;
        sortmode = 'ascend';
        if col(1)=='+'
            col = col(2:end);
        elseif col(1)=='-'
            sortmode = 'descend';
            col = col(2:end);
        end
        
        [~,idx] = sort(cellfun(@(C)mean(C),{result.(col)}),sortmode);
        result = result(idx);
    end
elseif nargin > 3
    % do the measurement
    
    nThread = getenv('OMP_NUM_THREADS');
    nThread = str2double(nThread);
    if isempty(nThread) || isnan(nThread)
        nThread = -1;
    end
    
    % check parallel pool
    pPool = gcp('nocreate');
    if isempty(pPool)
        nWorker = 1;
    else
        nWorker = pPool.NumWorkers;
    end
    
    if iscell(argin{2})
        argin{2} = sw_qscan(argin{2});
    end
    
    nQ = size(argin{2},2);
    if nQ == 3
        nQ = argin{2}(3);
    end
    
    if nargin > 3
        argin = [argin(:)' {'optmem' nMemSlice}];
    end
    if nargin<5
        nRun = 1;
    end
    % number of spin wave modes
    nMode = argin{1}.nmagext*2;
    % number of Q points triples for incommensurate structure
    if any(mod(argin{1}.mag_str.k(:,1),1))
        nQ = nQ*3;
    end
    
    swpref.setpref('usemex',usemex,'fid',0,'tid',0);
    
    tMeas = zeros(1,nRun);
    
    % measurement...
    for ii = 1:nRun
        tStart    = tic;
        spec      = fun(argin{:});
        tMeas(ii) = toc(tStart);
    end
    
    pps   = nQ./tMeas;
    
    if nRun == 1
        tMeas = round(tMeas*1e3)/1e3;
        pps   = round(pps*1e3)/1e3;
    end
    nMemSlice = spec.param.nSlice;
    hermit    = spec.param.hermit;
    
    result.time     = tMeas;
    result.fun      = fun;
    %result.argin    = argin;
    %result.argout   = spec;
    result.usemex   = usemex;
    result.nQ       = nQ;
    result.nWorker  = nWorker;
    result.nThread  = nThread;
    result.nSlice   = nMemSlice;
    result.pps      = pps;
    result.nRun     = nRun;
    result.nMode    = nMode;
    result.hermit   = hermit;
    % save result
    if exist(fName,'file')
        temp  = load(fName);
        result = [temp.result(:)' result];
        save(fName,'result');
    else
        save(fName,'result');
    end
    result = result(end);
end

for ii = 1:numel(result)
    rs       = result(ii);
    str01   = {'' 'x'};
    tMeasStr = err2str(mean(rs.time),std(rs.time));
    ppsStr   = err2str(mean(rs.pps),std(rs.pps));
    
    if rs.nThread<1
        nThreadStr = 'max';
    else
        nThreadStr = num2str(rs.nThread);
    end
    
    fprintf('%-30s | %-15s | %-10s | %-7d | %-7d | %-7d | %-7s | %-7d | %-7d | %-7s | %-7s\n',func2str(rs.fun),ppsStr,tMeasStr,...
        rs.nQ,rs.nMode,rs.nWorker,nThreadStr,rs.nSlice,rs.nRun,str01{rs.usemex+1},str01{rs.hermit+1});
end

end