function BiFeO3l()
    % setup
    swpref.setpref('usemex', true)
    bfo = spinw;
    bfo.genlattice('lat_const', [5.58 5.58 13.86], 'angled', [90 90 120], ...
                   'sym', 'R 3 c');
    bfo.addatom('r', [0 0 0.2212], 'S', 2.5, 'label', 'MFe3');
    bfo.gencoupling('maxDistance', 20)
    bfo.addmatrix('label', 'J1', 'value', 0.6904);
    bfo.addmatrix('label', 'J2', 'value', 0.2)
    bfo.addcoupling('mat', 'J1', 'bond', 1)
    bfo.addcoupling('mat', 'J2', 'bond', 2)
    bfo.addmatrix('label', 'D', 'value', [1 -1 0]*0.185)
    bfo.addcoupling('mat', 'D', 'bond', 2)
    bfo.optmagk('kbase', [1; 1; 0], 'seed', 1);
    
    % test parameters
    spinwave_args = {{[-1/2 0 0], [0 0 0], [1/2 1/2 0], 30}, ...
                     'sortMode', false, 'hermit', 1};
    egrid_args = {'component','Sperp','Evect',0:0.1:5};
    inst_args = {'dE',0.1};
    
    % do a spin wave calculation for incommensurate case ([nExt=[1,1,1]) 
    % and commensurate nExt=0.01 - which corresponds to a supercell of 
    % [11 11 1]
    for nExt = {[1, 1, 1], 0.01}
        bfo.genmagstr('nExt', nExt{1})
        for hermit = 0:1
            spinwave_args{end} = hermit;
            test_name = ['bifeo3_nExt_', regexprep(num2str(nExt{1}), ...
                                                   ' +', '_'), ...
                         '_hermit_', num2str(hermit)];
            sw_tests.utilities.profile_spinwave(test_name, bfo, ...
                                                spinwave_args, egrid_args, ...
                                                inst_args);
        end
    end
    


end