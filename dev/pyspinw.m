function pyspinw(varargin)
% the main application of the deployed SpinW app

fprintf('New SpinW process is created...\n')

switch nargin
    case 2
        % run transplant
        input = regexp(varargin{2},'(msgpack|json)'',''([\w:/]+)','tokens');
        msgformat = input{1}{1};
        url       = input{1}{2};
        
        transplant_remote(msgformat,url);
    case 1
        % execute command
        evalin('base',varargin{1});
end

end