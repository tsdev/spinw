function sw_notify(varargin)
% sends notification in OSX
%
% SW_NOTIFY(message)
%
% SW_NOTIFY('option1',value1, ...)
%
%
% The function can be used to show a notification when a calculation is
% finished. Note that the notification only shows up if Matlab is not
% the active window.
%
% Install terminal-notifier:
% Using Ruby, you can install it through RubyGems:
% $ [sudo] gem install terminal-notifier
% You can also install it via Homebrew:
% $ brew install terminal-notifier
%
% Options:
%


if ~ismac
    warning('sw_notify:WrongOS','The function only works in OSX!')
    return
end

path0 = '';

if nargin == 1
    param.message  = varargin{1};
    param.path     = path0;
    param.group    = 'spinw';
    param.subtitle = 'SpinW';
    param.sound    = 'ping';
else
    inpForm.fname  = {'message' 'path' 'group' 'subtitle' 'sound'};
    inpForm.defval = {'Done!'   path0  'spinw' 'SpinW'    'ping' };
    inpForm.size   = {[1 -1]    [1 -2]  [1 -3] [1 -4]     [1 -5] };
    inpForm.soft   = {false     true    false  false      false  };
    
    param = sw_readparam(inpForm, varargin{:});
end

cmd = [param.path 'terminal-notifier -message "' param.message '" -group "' param.group...
    '" -title "MATLAB" -sender com.mathworks.matlab -subtitle ' param.subtitle...
    ' -sound ' param.sound];

% execute bash_profile
[~,~] = system(['source ~/.bash_profile && ' cmd]);

end