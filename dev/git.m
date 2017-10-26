function git(varargin)
% basic git commands for SpinW repository
%
% ### Syntax
%
% `git('commit',msg)`
%
% `git('pull')`
%
% `branch = git('branch')`
%
% `git('branch',branch)`
%
% ### Description
%
% `git('commit',msg)` add all changes and commit with the given message.
%
% `git('pull')` pull updates from origin.
%
% `branch = git('branch')` show the current branch name.
%
% `git('branch',branch)` switch to `branch`.
%

if nargin > 0
    cmd = varargin{1};
end

if nargin>1
    msg = cellfun(@(C)[C ' '],varargin(2:end),'uniformoutput',false);
    msg = [msg{:}];
    msg = msg(1:(end-1));
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
    case 'pull'
        gitFun('pull');
    case 'branch'
        if nargin<2
            %[~,varargout{1}] = gitFun('rev-parse --abbrev-ref HEAD');
            gitFun('rev-parse --abbrev-ref HEAD');
        else
            branch = msg;
            gitFun('pull');
            gitFun(['checkout ' branch]);
        end
end

end