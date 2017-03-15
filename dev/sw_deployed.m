function sw_deployed(varargin)
% the main application of the deployed SpinW app

fprintf('Creating new SpinW process...\n')
sw_version;
input = regexp(varargin{2},'(msgpack|json)'',''([\w:/]+)','tokens');
msgformat = input{1}{1};
url       = input{1}{2};

transplant_remote(msgformat,url);

end