classdef jobqueue < handle
    % class to handle job submission server controlled over TCP/IP
    %    
    
    properties (Access = private)
        buffer  = struct('command',{},'jobID',{},'maxTime',{});
        % buffer structure that stores job information
        t       = [];
        % handle to the TCP/IP connection
        active = struct('command',{},'jobID',{},'maxTime',{});
        % stores the active job
    end
    
    methods (Static)
        function log(varargin)
            % write to STDOUT the given message with datetime
            fprintf([char(datetime) ' ' varargin{1} '\n'],varargin{2:end});
        end
    end
    
    methods
        
        function this = jobqueue(port)
            % create jobqueue by establishing TCP/IP connection to a client
            
            if port~=round(port) || port>65535 || port<1
                error('jobqueue:WrongInput','Port number has to be integer between 1-65535!')
            end
            
            this.log('TCP/IP connection initiated on port #%d',round(port));
            this.log('TCP/IP connection waiting for client to connect...');
            this.t = tcpip('0.0.0.0',port,'NetworkRole','server');
            fopen(this.t);
            this.log('TCP/IP channel open');
        end
        
        function str = readTCPIP(this)
            %
            % jobqueue.read()
            %
            % read string from TCPIP remove newline characters and convert
            % to horizontal character array
            
            str = char(fread(this.t, this.t.BytesAvailable)); %#ok<FREAD>
            str = str(str~=newline);
            str = str(:)';
            
        end
        
        function nJob = njob(this)
            % number of jobs pending
            
            nJob = numel(this.buffer);
        end
        
        function job = pop(this)
            % return the first job in the queue and register it as active
            
            if ~isempty(this.buffer)
                job = this.buffer(1);
                this.active = job;
                this.buffer = this.buffer(2:end);
            else
                job = this.buffer;
                this.active  = this.buffer;
            end
        end
        
        function done(this)
            % remove the active job
            
            this.active = this.active([]);
        end
        
        function signal = receive(this)
            % receive new job commands and add them into the queue
            %
            % signal = jobqueue.receive
            %
            % signal is a scalar with one of the following values:
            %
            % 0     No signal received.
            % 1     Some received (anything that is not 2)
            % 2     STOP signal received for the current job
            %
            
            signal = 0;
            if this.t.BytesAvailable>0
                signal = 1;
                % read the buffer
                tStart = tic;
                cmdTemp = this.readTCPIP;
                % read until the end command symbol reached ':'
                rTime = toc(tStart);
                while cmdTemp(end)~=':' && rTime<1
                    if this.t.BytesAvailable
                        cmdTemp = [cmdTemp this.readTCPIP]; %#ok<AGROW>
                    end
                    rTime = toc;
                end
                
                % split the commands
                cmdTemp = strsplit(cmdTemp,':');
                if rTime>=1
                    this.log('TCP/IP timed out (command end sign '':'' missing), throwing away last partial command ''%s''',cmdTemp{end});
                end
                
                % the last element in the cell is either a partial command or empty
                cmdTemp = cmdTemp(1:(end-1));
                if ~isempty(cmdTemp)
                    % keep only the valid EXE/STO commands
                    isExec = cellfun(@(C)~isempty(C),regexp(cmdTemp,'^EXEC \w+ [.0-9]+$'));
                    isStop = cellfun(@(C)~isempty(C),regexp(cmdTemp,'^STOP \w+$'));
                    isExit = cellfun(@(C)~isempty(C),regexp(cmdTemp,'^EXIT$'));
                    isValid = isExec|isStop|isExit;
                    
                    % list thrown away commands
                    if any(~isValid)
                        nisValid = find(~isValid);
                        for ii = 1:numel(nisValid)
                            this.log('Throwing away invalid command ''%s:''',cmdTemp{nisValid(ii)});
                        end
                    end
                    cmdTemp = cmdTemp(isValid);
                    
                    % convert jobid, etc
                    for ii = 1:numel(cmdTemp)
                        cmdTemp1 = strsplit(cmdTemp{ii},' ');
                        cmd   = cmdTemp1{1};
                        
                        switch cmd
                            case 'EXEC'
                                % only add the jobId if it is not identical to any
                                % existing job
                                jobID = cmdTemp1{2};
                                if ~ismember(jobID,{this.buffer.jobID})
                                    this.log('Added job ''%s'' to the execution queue',jobID);
                                    this.buffer(end+1).jobID = jobID;
                                    this.buffer(end).command = 'EXEC';
                                    this.buffer(end).maxTime = str2double(cmdTemp1{3});
                                else
                                    this.log('Job ''%s'' is already in the execution queue',jobID);
                                end
                            case 'STOP'
                                jobID = cmdTemp1{2};
                                % remove queued EXE commands that are stopped
                                if ~isempty(this.active) && jobID == this.active.jobID
                                    % stop signal received for the
                                    % currently running job, don't remove
                                    % it from the list
                                    signal = 2;
                                else
                                    if ismember(jobID,{this.buffer.jobID})
                                        this.log('Stopped job ''%s'' before execution',jobID);
                                    else
                                        this.log('Job ''%s'' is done before stop command sent',jobID);
                                    end
                                    this.buffer = this.buffer(~ismember({this.buffer.jobID},jobID));
                                end
                            case 'EXIT'
                                % register exit command
                                this.log('Exit signal registered');
                                this.buffer(end+1).command = 'EXIT';
                        end
                    end
                end
            end
        end
    end
end