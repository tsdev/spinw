function git(varargin)
% basic git commands
%
% ### Syntax
%
% `git commit Commit message`
%
% `git branch`
%
% `git branch branchName`
%
% `git command str1 str2 ...`
%
% ### Description
%
% `git commit Commit message` adds all changes and commit with the given
% message.
%
% `git pull` pull updates from origin.
%
% `git branch` show the current branch name.
%
% `git branch branchName` switch to `branch`.
%
% `git command str1 str2 ...` execute any git command with arbitrary number
% of strings to follow
%
% ### Examples
%
% Any git command can be used, for example:
%
% ```
% git add .
% git pull
% ```
%

if nargin > 0
    cmd = varargin{1};
end

if nargin>1
    msg = cellfun(@(C)[C ' '],varargin(2:end),'uniformoutput',false);
    msg = [msg{:}];
    msg = msg(1:(end-1));
else
    msg = '';
end

dir0   = pwd;
c      = onCleanup(@(~)cd(dir0));
cd(sw_rootdir);
gitFun = @(str)system(['git ' str],'-echo');

switch cmd
    case 'commit'
        if nargin<2
            msg = input('Commit message:');
        end
        if isempty(msg)
            warning('git:MissingCommitMessage','Commit message is required!');
        else
            gitFun('add .');
            gitFun(['commit -m "' msg '"']);
            gitFun('push');
        end
    case 'branch'
        if isempty(msg)
            %[~,varargout{1}] = gitFun('rev-parse --abbrev-ref HEAD');
            gitFun('rev-parse --abbrev-ref HEAD');
        else
            branch = msg;
            gitFun('pull');
            gitFun(['checkout ' branch]);
        end
    otherwise
        % any git command
        gitFun([cmd ' ' msg]);
end

end