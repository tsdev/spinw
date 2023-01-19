function profile_spinwave(test_name, sw_obj, spinwave_args, egrid_args, ...
                          inst_args, do_profile)

    if nargin < 6
        do_profile = true;
        % start profiling
        profile('clear');
        profile('on', '-memory');
    end
    % use supercell k=0 structure
    start_time = tic;
    spec = sw_obj.spinwave(spinwave_args{:});
    fprintf('Elapsed time for spinwave = %.4e seconds\n', toc(start_time))
    if ~isempty(egrid_args)
        tic;
        spec = sw_egrid(spec, egrid_args{:});
        fprintf('Elapsed time for sw_egrid = %.4e seconds\n', toc)
        if ~isempty(inst_args)
            tic;
            sw_instrument(spec, inst_args{:});
            fprintf('Elapsed time for sw_instrument = %.4e seconds\n', toc)
        end
    end
    fprintf('Total time elapsed = %.4e seconds\n', toc(start_time))
   
    if do_profile
        % save profile results
        p = profile('info');

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
        profsave(p, save_dir);  % will mkdir if not exist
        sw_tests.utilities.save_profile_results_to_txt(p, save_dir)
        profile('off');
    end
end