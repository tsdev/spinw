function run_performance_tests(commit)
    disp(version);
    if ~exist('spinw', 'class')
        if exist('swfiles', 'dir') && exist('external', 'dir') && exist('dat_files', 'dir')
            addpath(genpath('swfiles'));
            addpath(genpath('external'));
            addpath(genpath('dat_files'));
        else
            error(['SpinW is not installed and the swfiles, external and/or ', ...
                   'dat_files drectories couldn''t be found on the current ', ...
                   'path, so the tests cannot be run.'])
        end
    end
    
    % checkout a hard-coded commit sha for last updated benchmarking result
    % unless a specific commit (or branch name) is provided.
    if nargin == 0
        commit = "b0e60e32aa3790416b65386a721682e5361c1fb3";
    end
    current_branch = evalc("!git rev-parse --abbrev-ref HEAD");
    try
        eval(sprintf("!git checkout %s", commit));
    catch
        error('Could not checkout the commit provided.')
    end
    
    % run tests
    fpath_parts = {'+sw_tests', '+performance_tests'};
    files  = dir(fullfile(fpath_parts{:}, '*.m'));
    for fname = {files.name}
        run(erase(strjoin([fpath_parts, fname{1}], '.'), '+'));
    end

    % checkout initial branch
    eval(sprintf("!git checkout %s", current_branch(1:end-1)))
end
