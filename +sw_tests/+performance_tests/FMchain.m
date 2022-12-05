function FMchain()
    % setup
    swpref.setpref('usemex', true)
    FMchain = spinw;
    FMchain.genlattice('lat_const',[3 8 8],'angled',[90 90 90])
    FMchain.addatom('r', [0 0 0],'S', 1,'label','MCu1')
    FMchain.gencoupling('maxDistance',7)
    FMchain.addmatrix('value',-eye(3),'label','Ja')
    FMchain.addcoupling('mat','Ja','bond',1);
    FMchain.genmagstr('mode','direct', 'k',[0 0 0], ...
                      'n',[1 0 0],'S',[0; 1; 0]);

    
    % test parameters
    spinwave_args_common = {{[0 0 0], [1 0 0], 1e7}, ...
                            'sortMode', false};
    egrid_args = {'component','Sperp','Evect',0:0.1:5};
    inst_args = {'dE',0.1};
    

    for hermit = 0:1
        for optmem = 0:5:10
            spinwave_args = [spinwave_args_common, ...
                             'hermit', hermit, 'optmem', optmem];
            test_name = [mfilename() '_hermit_', num2str(hermit), ...
                         '_optmem_' num2str(optmem)];
            sw_tests.utilities.profile_spinwave(test_name, FMchain, ...
                                                spinwave_args, egrid_args, ...
                                                inst_args);
        end
    end
    


end