function profile_spinwave(test_name, sw_obj, spinwave_args, egrid_args, ...
                          inst_args)

    % start profiling
    profile('clear');
    profile('on', '-memory');
            
    % use supercell k=0 structure
    spec = sw_obj.spinwave(spinwave_args{:});
    if ~isempty(egrid_args)
        spec = sw_egrid(spec, egrid_args{:});
        if ~isempty(inst_args)
            sw_instrument(spec, inst_args{:});
        end
    end
   
    % save profile results
    p = profile('info');
    
    ver = sw_version();
    host_info = [computer(), '_', version('-release'), '_', ...
        ver.Name, ver.Release];
    save_dir = fullfile(pwd, "profile_results", host_info, test_name);
    profsave(p, save_dir);  % will mkdir if not exist
    sw_tests.utilities.save_profile_results_to_txt(p, save_dir)
    profile('off');
end