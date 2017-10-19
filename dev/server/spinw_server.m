function spinw_server(varargin)
% runs the SpinW server
%
% spinw_server path numWorkers portNum
%

fprintfd('SpinW Server started.\n');

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

if ~ischar(varargin{1}) || ~isnumeric(varargin{2}) || numel(varargin{2})~=1 ||...
        ~isnumeric(varargin{3}) || numel(varargin{3})~=1
    error('spinw_server:WrongInput',['Wrong input arguments!\nCall "spinw_server'...
        ' path numWorker portNum", e.g. "spinw_server /home/user1/srv 4 4002"!']);
end

folder  = varargin{1};
nWorker = varargin{2};
portNum = varargin{3};

% path to the log file
logPath = [folder filesep 'spinw_server_log.txt'];

if folder(end) == filesep
    % remove trailing slash
    folder = folder(1:(end-1));
end

fprintfd('TCP/IP connection initiated...\n');
% start the TCP/IP server to listen on the given port
t = tcpip('0.0.0.0',portNum,'NetworkRole','server');
fopen(t);
fprintfd('TCP/IP channel open.\n');

% start the parallel pool
if nWorker > 0
    fprintfd('Parallel pool starting up with %d worker...\n',nWorker);
    parpool(nWorker);
    fprintfd('Parallel pool started.\n');
else
    fprintfd('No parallel pool initiated, jobs will run on single process.\n');
end

buffer = struct('command',{},'jobID',{},'maxTime',{});

% start execution loop
while 1
    if t.BytesAvailable>0
        % read the buffer
        tic;
        cmdTemp = readTCPIP(t);
        % read until the end command symbol reached ':'
        rTime = toc;
        while cmdTemp(end)~=':' && rTime<1
            if t.BytesAvailable
                cmdTemp = [cmdTemp readTCPIP(t)]; %#ok<AGROW>
            end
            rTime = toc;
        end
        
        % split the commands
        cmdTemp = strsplit(cmdTemp,':');
        if rTime>=1
            fprintfd('TCP/IP timed out (command end sign ":" missing), throwing away last partial command "%s"!\n',cmdTemp{end});
        end
        
        % the last element in the cell is either a partial command or empty
        cmdTemp = cmdTemp(1:(end-1));
        if ~isempty(cmdTemp)
            % keep only the valid EXE/STO commands
            isExec  = cellfun(@(C)~isempty(C),regexp(cmdTemp,'^EXEC \w+ [.0-9]+$'));
            isStop  = cellfun(@(C)~isempty(C),regexp(cmdTemp,'^STOP \w+$'));
            isValid = isExec|isStop;
            
            % list thrown away commands
            if any(~isValid)
                nisValid = find(~isValid);
                for ii = 1:numel(nisValid)
                    fprintfd('Throwing away invalid command "%s:"!\n',cmdTemp{nisValid(ii)});
                end
            end
            cmdTemp = cmdTemp(isValid);
            
            % convert jobid, etc
            for ii = 1:numel(cmdTemp)
                cmdTemp1 = strsplit(cmdTemp{ii},' ');
                cmd   = cmdTemp1{1};
                jobID = cmdTemp1{2};
                
                switch cmd
                    case 'EXEC'
                        % only add the jobId if it is not identical to any
                        % existing job
                        if ~ismember(jobID,{buffer.jobID})
                            fprintfd('Added job "%s" to the execution queue!\n',jobID);
                            buffer(end+1).jobID = jobID; %#ok<AGROW>
                            buffer(end).command = 'EXE';
                            buffer(end).maxTime = str2double(cmdTemp1{3});
                        else
                            fprintfd('Job "%s" is already in the execution queue!\n',jobID);
                        end
                    case 'STOP'
                        % remove queued EXE commands that are stopped
                        if ismember(jobID,{buffer.jobID})
                            fprintfd('Stopped job "%s" before execution!\n',jobID);
                        end
                        buffer = buffer(~ismember({buffer.jobID},jobID));
                end
            end
        end
    end
    
    % excute the first command if the buffer is not empty
    if ~isempty(buffer)
        
        if ~exist([folder filesep 'in_' jobID '.mat'],'file')
            % wrong jobID move on
            fprintfd('Missing job "%s" source file!\n',buffer(1).jobID);
        elseif exist([folder filesep 'out_' jobID '.mat'],'file')
            fprintfd('Job "%s" output file already exists!\n',buffer(1).jobID);
        else
            % try to excute the command
            try
                fprintfd(['Executing job "%s" ...\n' repmat('-',1,80) '\n'],buffer(1).jobID);
                diary(logPath);
                diary('on');
                % load the .mat file
                input = load([folder filesep 'in_' buffer(1).jobID '.mat']);
                % excute the command
                if all(isfield(input,{'fun' 'argin' 'nargout'})) && ischar(input.fun) &&...
                        iscell(input.argin) && isnumeric(input.nargout) && numel(input.nargout)==1
                    if round(input.nargout)~=input.nargout || input.nargout<1
                        error('spinw_server:WrongNargOut','The nargout variable has to be integer and larger than zero!');
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
                    switch input.fun
                        case 'spinwave'
                            argout = cell(1,input.nargout);
                            [argout{:}] = spinwave(input.argin{:}); %#ok<NASGU>
                        case 'powspec'
                            argout = cell(1,input.nargout);
                            [argout{:}] = powspec(input.argin{:}); %#ok<NASGU>
                        otherwise
                            error('spinw_server:WrongFunction','The given function ''%s'' is not supported!',input.fun);
                    end
                else
                    error('spinw_server:WrongMatFile','The variables saved in the given .mat file has wrong format!');
                end
                diary('off');
                fprintf([repmat('-',1,80) '\n']);
                fprintfd('Job "%s" finished succesfully!\n',buffer(1).jobID);
                err = []; %#ok<NASGU>
            catch err %#ok<NASGU>
                argout = {}; %#ok<NASGU>
                diary('off')
                fprintf([repmat('-',1,80) '\n']);
                fprintfd('Job "%s" failed!\n',buffer(1).jobID);
            end
            
            % try to read the diary
            try
                log    = fileread(logPath); %#ok<NASGU>
            catch
                log = ''; %#ok<NASGU>
            end
            
            % delete the log file
            if exist(logPath,'file')
                delete(logPath);
            end
            % save the empty result with exception
            save([folder filesep 'out_' buffer(1).jobID '.mat'],'argout','err','log')
            % remove variables from memory
            clear('input','argout','err','log');
        end
        % remove the entry in the buffer
        buffer = buffer(2:end);

    end
end

end

function fprintfd(varargin)
% fprintf with date

fprintf([char(datetime) ' ' varargin{1}],varargin{2:end});

end

function str = readTCPIP(t)
% read string and remove newline character

str = char(fread(t, t.BytesAvailable)); %#ok<FREAD>
str = str(str~=newline);
str = str(:)';

end