classdef unittest_spinw_spinwave < sw_tests.unit_tests.unittest_super
    % Runs through unit test for @spinw/spinwave.m

    properties
        swobj = [];
        swobj_tri = [];
        default_spinwave = struct('formfact', false, ...
                                  'incomm', false, ...
                                  'helical', false, ...
                                  'norm', false, ...
                                  'nformula', int32(0), ...
                                  'param', struct('notwin', true, ...
                                                  'sortMode', true, ...
                                                  'tol', 1e-4, ...
                                                  'omega_tol', 1e-5, ...
                                                  'hermit', true), ...
                                  'title', 'Numerical LSWT spectrum', ...
                                  'gtensor', false, ...
                                  'datestart', '', ...
                                  'dateend', '');
        qh5 = [0:0.25:1; zeros(2,5)];
    end

    properties (TestParameter)
        % Test directions and literal qpts work
        qpts_h5 = {{[0 0 0], [1 0 0], 5}, ...
                   [0:0.25:1; zeros(2,5)]};
        mex = {0, 1};
    end

    methods (TestClassSetup)
        function setup_chain_model(testCase)
            % Just create a very simple FM 1D chain model
            testCase.swobj = spinw;
            testCase.swobj.genlattice('lat_const', [3 8 8], 'angled', [90 90 90]);
            testCase.swobj.addatom('r', [0 0 0],'S', 1, 'label', 'MNi2');
            testCase.swobj.gencoupling('maxDistance', 7);
            testCase.swobj.addmatrix('value', -eye(3), 'label', 'Ja');
            testCase.swobj.addcoupling('mat', 'Ja', 'bond', 1);
            testCase.swobj.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0; 1; 0]);
        end
        function setup_tri_model(testCase)
            % Create a simple triangular lattice model
            testCase.swobj_tri = spinw;
            testCase.swobj_tri.genlattice('lat_const', [4 4 6], 'angled', [90 90 120]);
            testCase.swobj_tri.addatom('r', [0 0 0], 'S', 3/2, 'label', 'MCr3');
            testCase.swobj_tri.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                                         'n', [0 0 1], 'k', [1/3 1/3 0]);
            J1 = 1;
            testCase.swobj_tri.addmatrix('label','J1','value',J1);
            testCase.swobj_tri.gencoupling;
            testCase.swobj_tri.addcoupling('mat','J1','bond',1);
        end
    end
    methods (TestMethodSetup)
        function disable_mex_setup(testCase)
            swpref.setpref('usemex', false);
        end
    end
    methods (TestMethodTeardown)
        function disable_mex_teardown(testCase)
            swpref.setpref('usemex', false);
        end
    end

    methods
        function sw = get_expected_sw_qh5(testCase)
            % Expected output for the chain model for 5 q-points from
            % [0 0 0] to [1 0 0]
            expected_hkl = testCase.qh5;
            expected_Sab = zeros(3, 3, 2, 5);
            Sab1 = [0.5 0  0.5j; 0 0 0; -0.5j 0 0.5];
            Sab2 = [0.5 0 -0.5j; 0 0 0;  0.5j 0 0.5];
            expected_Sab(:, :, 1, 1:4) = repmat(Sab1, 1, 1, 1, 4);
            expected_Sab(:, :, 2, 5) = Sab1;
            expected_Sab(:, :, 2, 1:4) = repmat(Sab2, 1, 1, 1, 4);
            expected_Sab(:, :, 1, 5) = Sab2;

            sw = testCase.default_spinwave;
            sw.omega = [ 1e-5  2.  4.  2. -1e-5; ...
                        -1e-5 -2. -4. -2.  1e-5];
            sw.Sab = expected_Sab;
            sw.hkl = expected_hkl;
            sw.hklA = expected_hkl*2/3*pi;
            sw.obj = copy(testCase.swobj);
        end
    end
    % Put tests with mocks in own block - prevent interference with other
    % tests
    methods (Test)
        function test_noInput(testCase)
            % Tests that if call spinwave with no input, it calls the help
            % First mock the help call
            help_function = sw_tests.utilities.mock_function('swhelp');
            testCase.swobj.spinwave();
            testCase.assertEqual(help_function.n_calls, 1);
            testCase.assertEqual(help_function.arguments, {{'spinw.spinwave'}});
        end
        function test_sw_qh5_optmem(testCase)
            qpts = testCase.qh5;
            optmem = 3;
            % Test that calculation is split into optmem chunks - sw_timeit
            % is called on each chunk plus once at beginning and end of
            % function. This is a bit fragile, may need to change the
            % target function, number of calls, or not check at all if
            % spinwave is refactored
            sw_timeit_mock = sw_tests.utilities.mock_function('sw_timeit');
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj.spinwave(qpts, 'optmem', optmem);
            testCase.assertEqual(sw_timeit_mock.n_calls, optmem + 2);
            % Test that with optmem gives the same result as without
            expected_sw = testCase.get_expected_sw_qh5();
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_sw_qh5_zero_freemem_warns(testCase)
            % Mock sw_freemem to return 0 to trigger warning
            sw_freemem_mock = sw_tests.utilities.mock_function( ...
                'sw_freemem', 0);
            sw_out = testCase.verifyWarning(...
                @() testCase.swobj.spinwave(testCase.qh5), ...
                'spinw:spinwave:FreeMemSize');
            testCase.verify_spinwave(sw_out, testCase.get_expected_sw_qh5);
        end
        function test_sw_qh5_low_freemem(testCase)
            % Check with low free memory calculation still attempted
            sw_freemem_mock = sw_tests.utilities.mock_function( ...
                'sw_freemem', 100);
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj.spinwave(testCase.qh5);
            testCase.verify_spinwave(sw_out, testCase.get_expected_sw_qh5);
        end
        function test_sw_qh5_fid(testCase)
            fprintf_mock = sw_tests.utilities.mock_function('fprintf0');
            fid = 3;
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj.spinwave(testCase.qh5, 'fid', fid);
            % check fid used to write file
            for irow = 1:fprintf_mock.n_calls
                testCase.assertEqual(fprintf_mock.arguments{irow}{1}, fid)
            end
            expected_sw = testCase.get_expected_sw_qh5();
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_sw_qh5_tid(testCase)
            sw_timeit_mock = sw_tests.utilities.mock_function('sw_timeit');
            tid = 2;
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj.spinwave(testCase.qh5, 'tid', tid);
            % check tid used in timing
            for irow = 1:sw_timeit_mock.n_calls
                testCase.assertEqual(sw_timeit_mock.arguments{irow}{3}, tid)
            end
            expected_sw = testCase.get_expected_sw_qh5();
            testCase.verify_spinwave(sw_out, expected_sw);
        end
    end
    methods (Test)
        function test_sw_qh5(testCase, qpts_h5, mex)
            swpref.setpref('usemex', mex);
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj.spinwave(qpts_h5);
            expected_sw = testCase.get_expected_sw_qh5();
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_sw_qh5_sortmode(testCase, mex)
            swpref.setpref('usemex', mex);
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj.spinwave(testCase.qh5, 'sortMode', false);
            expected_sw = testCase.get_expected_sw_qh5();
            % Sortmode swaps the last 2 modes
            expected_sw.omega([1 2], 5) = expected_sw.omega([2 1], 5);
            expected_sw.Sab(:, :, [1 2], 5) = expected_sw.Sab(:, :, [2 1], 5);
            expected_sw.param.sortMode = false;
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_sw_qh5_nformula(testCase, mex)
            swpref.setpref('usemex', mex);
            % Create copy to avoid changing obj for other tests
            swobj = copy(testCase.swobj);
            nformula = int32(2);
            swobj.unit.nformula = nformula;
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out_nformula = swobj.spinwave(testCase.qh5);
            expected_sw = testCase.swobj.spinwave(testCase.qh5);
            expected_sw.Sab = expected_sw.Sab/2;
            expected_sw.obj.unit.nformula = nformula;
            expected_sw.nformula = nformula;
            testCase.verify_spinwave(sw_out_nformula, expected_sw);
        end
        function test_sw_qh5_periodic(testCase, mex)
            swpref.setpref('usemex', mex);
            % Test qpts in different BZ give same omega, Sab
            qpts = testCase.qh5 + 1;
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj.spinwave(qpts);
            expected_sw = testCase.get_expected_sw_qh5();
            expected_sw.hkl = qpts;
            expected_sw.hklA =  [qpts(1, :)*2/3; qpts(2:end, :)*0.25 ]*pi;
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_sw_qh5_perpendicular(testCase, mex)
            swpref.setpref('usemex', mex);
            % Test qpts in perpendicular direction give flat modes
            qpts = [zeros(1, 5); 0:0.25:1; 0:0.25:1];
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj.spinwave(qpts);
            expected_sw = testCase.get_expected_sw_qh5();
            expected_sw.hkl = qpts;
            expected_sw.hklA =  [qpts(1, :)*2/3; qpts(2:end, :)*0.25 ]*pi;
            expected_sw.omega = 1e-5*[ones(1, 5); -ones(1, 5)];
            expected_sw.Sab(1, 3, :, 5) = -expected_sw.Sab(1, 3, :, 5);
            expected_sw.Sab(3, 1, :, 5) = -expected_sw.Sab(3, 1, :, 5);
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_sw_qh5_saveH_saveV(testCase, mex)
            swpref.setpref('usemex', mex);
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj.spinwave(testCase.qh5, ...
                                             'saveV', true, 'saveH', true);
            expected_V = repmat(eye(2), 1, 1, 5);
            expected_H = zeros(2, 2, 5);
            expected_H(:, :, [2 4]) = 2*repmat(eye(2), 1, 1, 2);
            expected_H(:, :, 3) = 4*eye(2);

            expected_sw = testCase.get_expected_sw_qh5();
            expected_sw.V = expected_V;
            expected_sw.H = expected_H;
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_sw_qh5_title(testCase, mex)
            swpref.setpref('usemex', mex);
            title = 'Example title';
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj.spinwave(testCase.qh5, 'title', title);
            expected_sw = testCase.get_expected_sw_qh5();
            expected_sw.title = title;
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_sw_with_nExt(testCase, mex)
            swpref.setpref('usemex', mex);
            % Create copy to avoid changing obj for other tests
            afm_chain = copy(testCase.swobj);
            afm_chain.matrix.mat = eye(3);
            afm_chain.genmagstr('mode', 'direct', 'k',[1/2 0 0], ...
                                'S',[0 0; 1 -1;0 0], 'nExt',[2 1 1]);
            sw_afm = afm_chain.spinwave(testCase.qh5);
            omega_vals = [0 2. 0 -2. 0];
            expected_omega = [omega_vals; omega_vals; -omega_vals; -omega_vals];
            testCase.verify_val(sw_afm.omega, expected_omega, 'abs_tol', 1e-7);
        end
        function test_sw_with_multiple_matom(testCase, mex)
            swpref.setpref('usemex', mex);
            fe_cu_chain = spinw;
            fe_cu_chain.genlattice('lat_const', [3 8 4], 'sym', 'P 1');
            fe_cu_chain.addatom('label', 'MCu2', 'r', [0 0 0]);
            fe_cu_chain.addatom('label', 'MFe2', 'r', [0 1/2 0]);
            fe_cu_chain.gencoupling;
            fe_cu_chain.addmatrix('label', 'J_{Cu-Cu}', 'value', 1);
            fe_cu_chain.addmatrix('label','J_{Fe-Fe}', 'value', 1);
            fe_cu_chain.addmatrix('label', 'J_{Cu-Fe}', 'value', -0.1)
            fe_cu_chain.addcoupling('mat','J_{Cu-Cu}','bond',1);
            fe_cu_chain.addcoupling('mat','J_{Fe-Fe}','bond',2);
            fe_cu_chain.addcoupling('mat','J_{Cu-Fe}','bond',[4 5]);
            fe_cu_chain.genmagstr('mode','helical','S',[0 0;1 1;0 0],'k',[1/2 0 0])

            testCase.disable_warnings('spinw:magstr:NotExact', 'spinw:spinwave:Twokm');
            sw_out = fe_cu_chain.spinwave(testCase.qh5, 'sortMode', false, 'omega_tol', 1e-12);
            om1 = 4.11473;
            om2 = 1.36015;
            om3 = 1.38527;
            expected_omega = zeros(12, 5);
            expected_omega(1:6, [1 3 5])= repmat([om2 0 0 -om2 om2 0]', 1, 3);
            expected_omega(1:6, [2 4])= repmat([om1 om3 -om3 -om1 om1 om3]', 1, 2);
            expected_omega(7:end, :) = -expected_omega(6:-1:1, :);
            testCase.verify_val(sw_out.omega, expected_omega, 'abs_tol', 5e-6);
        end
        function test_sw_saveSabp_commensurate_warns(testCase, mex)
            swpref.setpref('usemex', mex);
            sw = testCase.verifyWarning(...
                @() testCase.swobj.spinwave(testCase.qh5, 'saveSabp', true), ...
                'spinw:spinwave:CommensurateSabp');
            testCase.verify_spinwave(sw, testCase.get_expected_sw_qh5);
        end
        function test_sw_incom_in_supercell_warns(testCase, mex)
            swpref.setpref('usemex', mex);
            cycloid = spinw;
            cycloid.genlattice('lat_const',[3 8 10], 'sym',0);
            cycloid.addatom('r',[0 0 0],'S',1,'label','Cu1');
            cycloid.gencoupling('maxDistance',7);
            cycloid.addmatrix('label','J2','value', 1);
            cycloid.addcoupling('mat','J2','bond',2);
            % modulation of [1/4 0 0] gets transformed in genmagstr to
            % [0.5 0 0] for nExt = [2 1 1]
            testCase.disable_warnings('spinw:genmagstr:UCExtNonSuff');
            cycloid.genmagstr('mode', 'helical', ...
                              'S', [1 0; 0 1; 0 0], 'n', [0 0 1], ...
                               'nExt', [2 1 1], 'k', [0.25, 0, 0]);
            testCase.verifyWarning(@() cycloid.spinwave({[0 0 0], [1 0 0], 30}), ...
                                   {'spinw:spinwave:IncommKinSupercell', ...
                                   'spinw:spinwave:Twokm'});
        end
        function test_sw_qh5_saveSabp_incommensurate(testCase, mex)
            swpref.setpref('usemex', mex);
            qpts = testCase.qh5;
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj_tri.spinwave(qpts, ...
                                                 'saveSabp', true);
            expected_Sabp = zeros(3, 3, 2, 5);
            expected_Sabp(:, :, :, [1 5]) = repmat( ...
                diag([435.71079, 435.71079, 6.45497e-4]), 1, 1, 2, 2);
            expected_Sabp(:, :, :, [2 4]) = repmat( ...
                diag([0.59293, 0.59293, 0.47434]), 1, 1, 2, 2);
            expected_Sabp(:, :, :, 3) = repmat( ...
                diag([0.1875, 0.1875, 1.5]), 1, 1, 2, 1);
            omegap_vals =  [1.16190e-2 4.74342 3];
            expected_omegap = [ omegap_vals  omegap_vals(2:-1:1); ...
                               -omegap_vals -omegap_vals(2:-1:1)];

            testCase.verify_val(sw_out.Sabp, expected_Sabp, 'rel_tol', 1e-5);
            testCase.verify_val(sw_out.omegap, expected_omegap, 'rel_tol', 1e-5);
        end
        function test_sw_qh5_fitmode(testCase, mex)
            swpref.setpref('usemex', mex);
            qpts = testCase.qh5;
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj.spinwave(qpts, 'fitmode', true);
            % fitmode automatically turns off sortMode
            expected_sw = testCase.swobj.spinwave(qpts, 'sortMode', false);
            expected_sw = rmfield(expected_sw, {'obj', 'datestart', 'dateend'});
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_incommensurate(testCase, mex)
            swpref.setpref('usemex', mex);
            if mex
                err ='chol_omp:notposdef';
            else
                err= 'spinw:spinwave:NonPosDefHamiltonian';
            end
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            % Tests that incommensurate calculation is ok
            hkl = {[0 0 0] [0 1 0] [1 0 0] 5};
            % Create copy to avoid changing obj for other tests
            swobj = copy(testCase.swobj);
            commensurate_spec = swobj.spinwave(hkl);

            % Test incomm 2nd neighbour afm chain fails with only nearest
            % neighbour interactions
            swobj.genmagstr('mode', 'helical', 'k', [0.25 0 0], ...
                            'n', [0 0 1], 'S', [1; 0; 0]);
            testCase.verifyError(@() swobj.spinwave(hkl), err);
            % Add 2nd neighbour interactions and test this incommensurate
            % structure produces 3x number of modes as commensurate
            swobj.addmatrix('value', 2, 'label', 'Jb');
            swobj.addcoupling('mat', 'Jb', 'bond', 2);
            incomm_spec = swobj.spinwave(hkl);
            testCase.assertEqual(size(incomm_spec.omega, 1), ...
                                 size(commensurate_spec.omega, 1) * 3);
        end
        function test_twin(testCase, mex)
            swpref.setpref('usemex', mex);
            % Tests that setting twins gives correct outputs
            % Create copy to avoid changing obj for other tests
            swobj_twin = copy(testCase.swobj);
            swobj_twin.addtwin('axis', [0 0 1], 'phid', [60 120], 'vol', [1 2]);
            rotc = swobj_twin.twin.rotc;
            hkl = [1 2; 3 4; 5 6];
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = swobj_twin.spinwave(hkl);

            expected_sw = testCase.default_spinwave;
            expected_sw.param.notwin = false;
            expected_sw.omega = {};
            expected_sw.Sab = {};
            expected_sw.obj = copy(swobj_twin);
            expected_sw.hkl = hkl;
            expected_sw.hklA = [2/3 4/3; 0.75, 1; 1.25, 1.5]*pi;
            % Recalculate without twins for each set of hkl's and compare
            [qTwin, rotQ] = swobj_twin.twinq(hkl);
            for itwin = 1:3
                sw_single = testCase.swobj.spinwave(qTwin{itwin});
                expected_sw.omega = [expected_sw.omega sw_single.omega];
                rot = rotc(:, :, itwin);
                rot_Sab = zeros(3, 3, 2, 2);
                % Twin Sab is a rotation of single Sab
                for imode = 1:2
                    for iq = 1:2
                        rot_Sab(:, :, imode, iq) = rot*sw_single.Sab(:, :, imode, iq)*rot';
                    end
                end
                expected_sw.Sab = [expected_sw.Sab rot_Sab];
            end
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_notwin(testCase, mex)
            swpref.setpref('usemex', mex);
            % Create copy to avoid changing obj for other tests
            swobj_twin = copy(testCase.swobj);
            swobj_twin.addtwin('axis', [0 0 1], 'phid', [60 120], 'vol', [1 2]);
            qpts = testCase.qh5;
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = swobj_twin.spinwave(qpts, 'notwin', true);
            % Test even if twin is added to structure it is not actually
            % calculated if notwin is specified
            expected_sw = testCase.get_expected_sw_qh5;
            expected_sw.obj.twin = swobj_twin.twin;
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_cmplxBase_equivalent_with_tri(testCase, mex)
            swpref.setpref('usemex', mex);
            % For this structure, both cmplxBase true and false give the
            % same e-vectors so should give the same result
            qpts = {[0 0 0], [1 0 0], 5};
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj_tri.spinwave(qpts, 'cmplxBase', false);
            sw_out_cmplx = testCase.swobj_tri.spinwave(qpts, 'cmplxBase', true);
            testCase.verify_spinwave(sw_out_cmplx, sw_out);
        end
        function test_cmplxBase_fails_with_chain(testCase, mex)
            swpref.setpref('usemex', mex);
            if mex
                err ='chol_omp:notposdef';
            else
                err= 'spinw:spinwave:NonPosDefHamiltonian';
            end
            % Test cmplxBase actually does something - it should fail with
            % chain
            qpts = {[0 0 0], [1 0 0], 5};
            testCase.verifyError(...
                @() testCase.swobj.spinwave(qpts, 'cmplxBase', true), ...
                err);
        end
        function test_formfact(testCase, mex)
            swpref.setpref('usemex', mex);
            qpts = {[0 0 0] [10 5 1] 19};
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_ff = testCase.swobj.spinwave(qpts, 'formfact', true);
            % Test that Sab with the form factor (ff) is explicitly the
            % same as Sab with no ff multiplied by ff
            % ff calculated with sw_mff and the scaling is F(Q)^2.
            expected_sw = testCase.swobj.spinwave(qpts, 'formfact', false);
            ff = sw_mff(testCase.swobj.unit_cell.label{1}, sw_ff.hklA);
            expected_sw.Sab = expected_sw.Sab.*permute(ff.^2, [1 3 4 2]);
            expected_sw.formfact = true;

            testCase.verify_spinwave(sw_ff, expected_sw);
        end
        function test_formfactfun(testCase, mex)
            swpref.setpref('usemex', mex);
            function F = formfactfun(atom_label, Q)
                F = sum(Q, 1);
            end
            qpts = {[0 0 0] [10 5 1] 19};
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_ff = testCase.swobj.spinwave(qpts, 'formfact', true, ...
                                            'formfactfun', @formfactfun);
            % Test that Sab with the form factor (ff) is explicitly the
            % same as Sab with no ff multiplied by ff
            expected_sw = testCase.swobj.spinwave(qpts, 'formfact', false);
            ff = formfactfun(testCase.swobj.unit_cell.label{1}, sw_ff.hklA);
            expected_sw.Sab = expected_sw.Sab.*permute(ff.^2, [1 3 4 2]);
            expected_sw.formfact = true;

            testCase.verify_spinwave(sw_ff, expected_sw, 'rel_tol', 1e-15);
        end
        function test_gtensor(testCase, mex)
            swpref.setpref('usemex', mex);
            qpts = {[0 0 0], [1 1 1], 5};
            gmat = diag([1, 2, 3]);
            % Create copy to avoid changing obj for other tests
            swobj_g = copy(testCase.swobj);
            swobj_g.addmatrix('label','g_1','value', gmat)
            swobj_g.addg('g_1')
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_g = swobj_g.spinwave(qpts, 'gtensor', true);
            % Also check that it warns that gtensor is not being used
            expected_sw = testCase.verifyWarning(...
                @() swobj_g.spinwave(qpts, 'gtensor', false), ...
                'spinw:spinwave:NonZerogTensor');
            expected_sw.Sab = expected_sw.Sab.*[1 2 3; 2 4 6; 3 6 9];
            expected_sw.gtensor = true;
            testCase.verify_spinwave(sw_g, expected_sw);
        end
        function test_gtensor_incomm(testCase, mex)
            swpref.setpref('usemex', mex);
            qpts = {[0 0 0], [1 1 1], 5};
            gmat = diag([1, 2, 3]);
            % Create copy to avoid changing obj for other tests
            swobj_g = copy(testCase.swobj_tri);
            swobj_g.addmatrix('label','g_1','value', gmat)
            swobj_g.addg('g_1')
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_g = swobj_g.spinwave(qpts, 'gtensor', true);
            % Check that Sab with g is same as Sab without g but multiplied
            % by g in each direction
            expected_sw = swobj_g.spinwave(qpts);
            expected_sw.Sab = expected_sw.Sab.*[2.25 2.25 4.5; ...
                                                2.25 2.25 4.5; ...
                                                4.5  4.5  9];
            expected_sw.gtensor = true;
            testCase.verify_spinwave(sw_g, expected_sw, 'rel_tol', 1e-15);
        end
        function test_hermit(testCase, mex)
            swpref.setpref('usemex', mex);
            % Create copy to avoid changing obj for other tests
            swobj_h = copy(testCase.swobj);
            % Tests that the 'hermit' option to switch to a non-hermitian calculation works
            % First make the model non-Hermitian by adding a large axial SIA perpendicular to the spins
            swobj_h.addmatrix('label', 'K', 'value', diag([-1 0 0]));
            swobj_h.addaniso('K');
            hkl = {[0 0 0] [0 1 0] [1 0 0] 50};
            % Check that calling it with 'hermit' on gives an error
            testCase.assertError(@() swobj_h.spinwave(hkl, 'hermit', true), ?MException);
            % Now check that there are imaginary eigenvalues/energies in the output with 'hermit' off
            spec = swobj_h.spinwave(hkl, 'hermit', false);
            testCase.assertGreaterThan(sum(abs(imag(spec.omega(:)))), 0);
        end
        function test_sw_qh5_tol(testCase, mex)
            swpref.setpref('usemex', mex);
            tol = 5e-4;
            qpts = testCase.qh5;
            swobj_tol = copy(testCase.swobj);
            % Generate magstr that deviates slightly from commensurate
            swobj_tol.genmagstr('mode', 'helical', 'k', [tol 0 0], ...
                                'n', [0 0 1], 'S', [0; 1; 0]);
            if ~mex
                % Check that without tol it is incommensurate - causes error
                testCase.verifyError(...
                    @() swobj_tol.spinwave(qpts), ...
                    'spinw:spinwave:NonPosDefHamiltonian');
            end
            % Check that with tol is approximated to commensurate
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = swobj_tol.spinwave(qpts, 'tol', tol);
            expected_sw = testCase.get_expected_sw_qh5;
            expected_sw.obj = swobj_tol;
            expected_sw.param.tol = tol;
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_sw_qh5_omega_tol(testCase, mex)
            swpref.setpref('usemex', mex);
            if mex
                err ='chol_omp:notposdef';
            else
                err= 'spinw:spinwave:NonPosDefHamiltonian';
            end
            qpts = testCase.qh5;
            % Check that with no added omega_tol Hamiltonian isn't positive
            % definite - causes error
            testCase.verifyError(...
                @() testCase.swobj.spinwave(qpts, 'omega_tol', 0), ...
                err);
            % Check with omega_tol the omega is changed appropriately
            omega_tol = 1;
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.swobj.spinwave(qpts, 'omega_tol', omega_tol);
            expected_sw = testCase.get_expected_sw_qh5;
            expected_sw.omega(:, 1) = [omega_tol -omega_tol];
            expected_sw.omega(:, end) = [-omega_tol omega_tol];
            expected_sw.param.omega_tol = omega_tol;
            testCase.verify_spinwave(sw_out, expected_sw);
        end
        function test_biquadratic_with_incomm_causes_error(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            % Set coupling to biquadratic
            swobj_tri.coupling.type(:) = 1;
            testCase.verifyError(...
                @() swobj_tri.spinwave(testCase.qh5), ...
                'spinw:spinwave:Biquadratic');
        end
        function test_no_magstr_causes_error(testCase)
            swobj = spinw;
            swobj.genlattice('lat_const', [3 8 8], 'angled', [90 90 90]);
            swobj.addatom('r', [0 0 0]);
            testCase.verifyError(...
                @() swobj.spinwave(testCase.qh5), ...
                'spinw:spinwave:NoMagneticStr');
        end
    end
    methods (Test, TestTags = {'Symbolic'})
        function test_sw_symbolic_no_qpts(testCase)
            swobj = copy(testCase.swobj);
            swobj.symbolic(true);
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            sw_out = testCase.verifyWarning(@() swobj.spinwave(), 'spinw:spinwave:MissingInput');

            symstr = '-Ja*exp(-pi*h*2i)*(exp(pi*h*2i) - 1)^2';
            expected_sw.ham = [str2sym(symstr) sym(0); ...
                               sym(0) str2sym(symstr)];
            expected_sw.omega = [str2sym(symstr(2:end)); str2sym(symstr)];
            expected_sw.obj = swobj;
            expected_sw.datestart = '';
            expected_sw.dateend = '';
            expected_sw.title = 'Symbolic LSWT spectrum';
            testCase.verify_spinwave(sw_out, expected_sw);
        end
    end
end
