classdef unittest_spinw_spec2MDHisto < sw_tests.unit_tests.unittest_super
    properties
    end
    methods (Test)
        function test_1_0_0(testCase)
            q0 = [0 0 0];
            qdir = [1 0 0];
            spec = sw_egrid(spinwave(sw_model('triAF', 1), {q0 q0+qdir nsteps}));
            proj = [qdir(:) [1 -1 0]' [0 0 1]'];
            dproj = [norm((qdir-q0))/nsteps, 1e-6, 1e-6];
            sw_spec2MDHisto(spec, proj, dproj, 'test100mdh.nxs');
        end

        function test_1_1_0(testCase)
            q0 = [0 0 0];
            qdir = [1 1 0];
            spec = sw_egrid(spinwave(sw_model('triAF', 1), {q0 q0+qdir nsteps}));
            proj = [qdir(:) [1 -1 0]' [0 0 1]'];
            dproj = [norm((qdir-q0))/nsteps, 1e-6, 1e-6];
            sw_spec2MDHisto(spec, proj, dproj, 'test110mdh.nxs');
        end

        function test_1_1_1(testCase)
            q0 = [0 0 0];
            qdir = [1 1 1];
            spec = sw_egrid(spinwave(sw_model('triAF', 1), {q0 q0+qdir nsteps}));
            proj = [qdir(:) [1 -1 0]' [1 1 -2]'];
            dproj = [norm((qdir-q0))/nsteps, 1e-6, 1e-6];
            sw_spec2MDHisto(spec, proj, dproj, 'test111mdh.nxs');
        end
        function test_1_1_2(testCase)
            q0 = [0 0 2];
            qdir = [1 1 0];
            spec = sw_egrid(spinwave(sw_model('triAF', 1), {q0 q0+qdir nsteps}));
            proj = [qdir(:) [1 -1 0]' [0 0 1]'];
            dproj = [norm((qdir-q0))/nsteps, 1e-6, 1e-6];
            sw_spec2MDHisto(spec, proj, dproj, 'tmp/test112mdh.nxs');
        end
    end