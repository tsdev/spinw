function git(varargin)
% basic git commands
%
% ### Syntax
%
% `git repo commit Commit message`
%
% `git repo branch`
%
% `git repo branch branchName`
%
% `git repo command str1 str2 ...`
%
% ### Description
%
% `git repo commit Commit message` adds all changes and commit with the
% given message to the given repo. The give repo is reached via the
% `go(repo)` command, thus repo can be the name of a Matlab function within
% the repository or a label that the `go` function knows about or a valid
% path.
%
% `git repo pull` pull updates from origin.
%
% `git repo branch` show the current branch name.
%
% `git repo branch branchName` switch to `branch`.
%
% `git repo command str1 str2 ...` execute any git command with arbitrary
% number of strings to follow
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
% ### See Also
%
% [go]
%

if nargin < 2
    return
end

package = varargin{1};
cmd     = varargin{2};

if nargin>2
    msg = cellfun(@(C)[C ' '],varargin(3:end),'uniformoutput',false);
    msg = [msg{:}];
    msg = msg(1:(end-1));
else
    msg = '';
end

dir0   = pwd;
c      = onCleanup(@(~)cd(dir0));
go(package);
gitFun = @(str)system(['git ' str],'-echo');

fprintf('Repository: ');
gitFun('rev-parse --show-toplevel');

switch cmd
    case 'commit'
        if nargin<2
            msg = input('Commit message:');
        end
        if isempty(msg)
            warning('git:MissingCommitMessage','Commit message is required!');
        else
            gitFun('add -A');
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