function profile_spinwave(test_name, sw_obj, spinwave_args, egrid_args, ...
                          inst_args, do_profile)

    % generate file directory name for results
    try
        % get current commit is possible
        commit = evalc("!git rev-parse --short HEAD");
        commit = commit(1:end-1); % remove newline
    catch
        ver = sw_version();
        commit = [ver.Name ver.Release];
    end
    host_info = [computer(), '_', version('-release')];
    save_dir = fullfile(pwd, "profile_results", commit, ...
        host_info, test_name);
    mkdir(save_dir)

    % open file for tic/toc timings
    fid = fopen(fullfile(save_dir, ...
                sprintf('tictoc_times_profile_%i.txt', do_profile)), 'w');
    c = onCleanup(@()fclose(fid));  % in case of exception file will close
    fprintf(fid, "Function\tDuration (s)\n");

    if nargin < 6
        do_profile = true;
        % start profiling
        profile('clear');
        profile('on', '-memory');
    end
    % use supercell k=0 structure
    start_time = tic;
    spec = sw_obj.spinwave(spinwave_args{:});
    fprintf(fid, 'spinwave\t%.4e\n', toc(start_time));
    if ~isempty(egrid_args)
        tic;
        spec = sw_egrid(spec, egrid_args{:});
        fprintf(fid, 'sw_egrid\t%.4e\n', toc);
        if ~isempty(inst_args)
            tic;
            sw_instrument(spec, inst_args{:});
            fprintf(fid, 'sw_instrument\t%.4e\n', toc);
        end
    end

    if do_profile
        % save profile results
        p = profile('info');
        profsave(p, save_dir);  % will mkdir if not exist
        % save ascii summary
        sw_tests.utilities.save_profile_results_to_txt(p, save_dir);
        profile('off');
    end
end