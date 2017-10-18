function spinw_server(varargin)
% runs the SpinW server
%
% spinw_server path numWorkers portNum
%

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

if ~ischar(varargin{1}) || ~isnumeric(varargin{2}) || numel(varargin{2})~=1 || ~isnumeric(varargin{3}) || numel(varargin{3})~=1
    error('spinw_server:WrongInput','Wrong input arguments!\nCall "spinw_server path numWorker portNum", e.g. "spinw_server /home/user1/srv 4 4002"!');
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

% start the TCP/IP server to listen on the given port
t = tcpip('0.0.0.0',portNum,'NetworkRole','server');
fopen(t);

% start the parallel pool
if nWorker > 0
    parpool(nWorker);
end

buffer = struct('command',{},'jobID',{},'maxTime',{});

% start execution loop
while 1
    if t.BytesAvailable>0
        % read the buffer
        cmdTemp = char(fread(t, t.BytesAvailable)); %#ok<FREAD>
        % read until the end command symbol reached ':'
        while cmdTemp(end)~=':'
            if t.BytesAvailable
                cmdTemp = [cmdTemp; char(fread(t, t.BytesAvailable))]; %#ok<AGROW,FREAD>
            end
        end
        
        % split the commands
        cmdTemp = strsplit(cmdTemp(:)',':');
        % keep only the valid EXE/STO commands
        isExe   = cellfun(@(C)~isempty(C),regexp(cmdTemp,'^EXE \w+ [.0-9]+$'));
        isStop  = cellfun(@(C)~isempty(C),regexp(cmdTemp,'^STO \w+$'));
        cmdTemp = cmdTemp(isExe|isStop);
        
        % convert jobid, etc
        for ii = 1:numel(cmdTemp)
            cmdTemp1 = strsplit(cmdTemp{ii},' ');
            cmd   = cmdTemp1{1};
            jobID = cmdTemp1{2};
            
            switch cmd
                case 'EXE'
                    % only add the jobId if it is not identical to any
                    % existing job
                    if ~ismember(jobID,{buffer.jobID})
                        buffer(end+1).jobID = jobID; %#ok<AGROW>
                        buffer(end).command = 'EXE';
                        buffer(end).maxTime = str2double(cmdTemp1{3});
                    end
                case 'STO'
                    % remove queued EXE commands that are stopped
                    buffer = buffer(~ismember({buffer.jobID},jobID));
            end
        end
    end
    
    % excute something is the buffer is not empty
    if ~isempty(buffer)
        % try to excute the command
        try
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
                    pPool = gcp('nocreate');
                    if isempty(pPool) && nWorker > 0
                        parpool(nWorker);
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
            err = []; %#ok<NASGU>
        catch err %#ok<NASGU>
            argout = {}; %#ok<NASGU>
        end
        
        % try to read the diary
        try
            log    = fileread(logPath); %#ok<NASGU>
        catch
            log = ''; %#ok<NASGU>
        end
        
        % delete the log file
        delete(logPath);
        % save the empty result with exception
        save([folder filesep 'out_' buffer(1).jobID '.mat'],'argout','err','log')
        % remove the entry in the buffer
        buffer = buffer(2:end);
        % remove variables from memory
        clear('input','argout','err','log');
    end
    
end

end







