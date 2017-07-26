function pyspinw(varargin)
% the main application of the deployed SpinW app

fprintf('SpinW process is running...\n')

switch nargin
    case 1
        % execute command
        evalin('base',varargin{1});

    otherwise
                
        % run transplant
        [~,input] = regexp(varargin{end},'transplant_remote\((.*?)\)','match','tokens');
        input     = strsplit(input{1}{1},',');
        %input = regexp(varargin{2},'(msgpack|json)'',''([\w:/]+)','tokens');
        msgformat = input{1}(2:end-1);
        url       = input{2}(2:end-1);
        %zmqpath   = input{3}(2:end-1);
        
        %transplant_remote(msgformat,url,zmqpath);
        switch computer
            case 'MACI64'
                transplant_remote(msgformat,url,'libzmq.dylib');
            case 'PCWIN64'
                transplant_remote(msgformat,url,'libzmq.dll');
            case 'GLNXA64'
                transplant_remote(msgformat,url,'libzmq.so');
            otherwise
                error('pyspinw:WrongArchitecture','Unsupported architecture!')
        end
        
end

end