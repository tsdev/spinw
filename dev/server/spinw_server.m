function spinw_server(varargin)
% runs the SpinW server
%
% spinw_server path numWorkers portNum
%

jobqueue.log('SpinW Server started0');

if nargin~=3
    fprintf('Input arguments:\n');
    for ii = 1:nargin
        argin = varargin{ii};
        if ischar(argin)
            fprintf('#%d: string ''%s''\n',ii,argin);
        elseif isnumeric(argin)
            fprintf('#%d: numeric [%dx%d]\n',ii,size(argin,1),size(argin,2));
        end
    end
    return
end

folder  = varargin{1};
nWorker = varargin{2};
nPort   = varargin{3};

% list of functions that can be executed
funList = {'spinwave' 'powspec' 'powspecfast' 'spinwavefast'};
funProf = [true       true      true          true          ];

if isdeployed
    % convert strings into double
    nWorker = str2double(nWorker);
    nPort = str2double(nPort);
end

if ~ischar(folder) || ~isnumeric(nWorker) || numel(nWorker)~=1 ||...
        ~isnumeric(nPort) || numel(nPort)~=1 || isnan(nWorker) || isnan(nPort)
    error('spinw_server:WrongInput',['Wrong input arguments!\nCall "spinw_server'...
        ' path numWorker portNum", e.g. "spinw_server /home/user1/srv 4 4002"!']);
end

% path to the log file
logPath = [folder filesep 'spinw_server_log.txt'];

if folder(end) == filesep
    % remove trailing slash
    folder = folder(1:(end-1));
end

queue = jobqueue(nPort);

% start the parallel pool
if nWorker > 0
    queue.log('Parallel pool starting up with %d worker...',nWorker);
    parpool(nWorker);
    queue.log('Parallel pool started');
else
    queue.log('No parallel pool initiated, jobs will run on single process');
end

isWaiting = false;
shutDown  = false;

% start execution loop
while ~shutDown
    if queue.receive > 0
        isWaiting = false;
    end
    
    % excute the first command if the buffer is not empty
    if queue.njob > 0
        isWaiting = false;
        
        % get the first job in the queue
        job = queue.pop;
        switch job.command
            case 'EXEC'
                fileIn  = [folder filesep 'in_' job.jobID '.mat'];
                fileOut = [folder filesep 'out_' job.jobID '.mat'];
                
                if ~exist(fileIn,'file')
                    % wrong jobID move on
                    queue.log('Missing job "%s" source file',job.jobID);
                elseif exist(fileOut,'file')
                    queue.log('Job "%s" output file already exists',job.jobID);
                else
                    % try to excute the command
                    try
                        queue.log(['Executing job "%s"...\n' repmat('-',1,80)],job.jobID);
                        diary(logPath);
                        diary('on');
                        % load the .mat file
                        try
                            input = load(fileIn);
                        catch
                            % OK, this is byte array, not a real mat file....
                            f = fopen(fileIn);
                            c = onCleanup(@() fclose(f));
                            input = getArrayFromByteStream(fread(f,Inf,'uint8=>uint8'));
                            queue.log('Reading input file as a byte array');
                        end
                        % excute the command
                        if all(isfield(input,{'fun' 'argin' 'nargout'})) && ischar(input.fun) &&...
                                iscell(input.argin) && isnumeric(input.nargout) && numel(input.nargout)==1
                            if round(input.nargout)~=input.nargout || input.nargout<1
                                error('spinw_server:WrongNargOut','The nargout variable has to be integer and larger than zero!');
                            end
                            
                            % get debug
                            if isfield(input,'prof') && isnumeric(input.prof) && ...
                                    numel(input.prof) == 1 && ismember(input.prof,[0 1 2])
                                prof.level = input.prof;
                                prof.argin = {};
                            elseif isfield(input,'prof') && iscell(input.prof) && isnumeric(input.prof{1}) && ...
                                    numel(input.prof{1}) == 1 && ismember(input.prof{1},[0 1 2])
                                prof.level = input.prof{1};
                                prof.argin = input.prof{2:end};
                            else
                                prof.level = 0;
                                prof.argin = {};
                            end
                            
                            % restart parallel pool if it timed out
                            if ismember(input.fun,{'spinwave' 'powspec'})
                                if nWorker > 0
                                    pPool = gcp('nocreate');
                                    if isempty(pPool)
                                        parpool(nWorker);
                                    end
                                end
                            end
                            
                            % run only the registered functions
                            funIdx = ismember(funList,input.fun);
                            if ~funProf(funIdx) && prof.level>1 && nWorker~=0
                                % functions that does not support mpiprofile
                                % switch to time measurement
                                prof.level = 1;
                            end
                            
                            if any(funIdx)
                                fun = str2func(funList{funIdx});
                                % call the function
                                argout      = cell(1,input.nargout);
                                switch prof.level
                                    case 0
                                        % no profging
                                        [argout{:}] = fun(input.argin{:}); %#ok<NASGU>
                                        prof.info  = [];
                                    case 1
                                        % measure execution time
                                        tStart      = tic;
                                        [argout{:}] = fun(input.argin{:}); %#ok<NASGU>
                                        prof.info.time = toc(tStart);
                                    case 2
                                        % try to do profiling
                                        if nWorker == 0
                                            % profile here
                                            profile('on',prof.argin{:});
                                            [argout{:}] = fun(input.argin{:}); %#ok<NASGU>
                                            prof.info = profile('info');
                                        else
                                            % mpiprofile within the funciton
                                            [argout{:},prof.info] = fun(input.argin{:},'prof',true); %#ok<ASGLU>
                                        end
                                end
                            else
                                error('spinw_server:WrongFunction','The given function ''%s'' is not supported!',input.fun);
                            end
                        else
                            error('spinw_server:WrongMatFile','The variables saved in the given .mat file has wrong format!');
                        end
                        diary('off');
                        fprintf([repmat('-',1,80) '\n']);
                        queue.log('Job "%s" finished succesfully',job.jobID);
                        err = []; %#ok<NASGU>
                    catch err %#ok<NASGU>
                        argout = {}; %#ok<NASGU>
                        prof = struct('level',{},'argin',{},'info',{});
                        diary('off')
                        fprintf([repmat('-',1,80) '\n']);
                        queue.log('Job "%s" failed',job.jobID);
                    end
                    
                    % try to read the diary and save it into the log
                    try
                        log    = fileread(logPath); %#ok<NASGU>
                    catch
                        log = ''; %#ok<NASGU>
                    end
                    
                    % delete the log file
                    if exist(logPath,'file')
                        delete(logPath);
                    end
                    % save the results, exception stack and log
                    if ~isempty(prof) && ~isempty(prof.info)
                        save(fileOut,'argout','err','log','prof');
                    else
                        save(fileOut,'argout','err','log');
                    end
                    % remove variables from memory
                    clear('input','argout','err','log');
                    % adding a '.done' at the end of the input file to signal its finished
                    movefile(fileIn,[fileIn '.done']);
                end
                
                % register that the active job is done
                queue.done;
            case 'EXIT'
                shutDown = true;
                queue.log('Server is shutting down...');
        end
    elseif ~isWaiting
        queue.log('Idle, waiting for new job...');
        isWaiting = true;
    end
end

end