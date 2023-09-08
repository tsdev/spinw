classdef unittest_spinw_spec2MDHisto < sw_tests.unit_tests.unittest_super
    
    properties(TestParameter)
        nsteps={100};
    end
    methods (Test)
        function test_1_0_0(testCase)
            q0 = [0 0 0];
            qdir = [1 0 0];
            spec = sw_egrid(spinwave(sw_model('triAF', 1), {q0 q0+qdir testCase.nsteps{1}}));
            proj = [qdir(:) [0 1 0]' [0 0 1]'];
            dproj = [norm((qdir-q0))/testCase.nsteps{1}, 1e-6, 1e-6];
            sw_spec2MDHisto(spec, proj, dproj, 'tmp/test100mdh.nxs');
        end

        function test_1_1_0(testCase)
            q0 = [0 0 0];
            qdir = [1 1 0];
            spec = sw_egrid(spinwave(sw_model('triAF', 1), {q0 q0+qdir testCase.nsteps{1}}));
            proj = [qdir(:) [1 -1 0]' [0 0 1]'];
            dproj = [norm((qdir-q0))/testCase.nsteps{1}, 1e-6, 1e-6];
            sw_spec2MDHisto(spec, proj, dproj, 'tmp/test110mdh.nxs');
        end

        function test_1_1_1(testCase)
            q0 = [0 0 0];
            qdir = [1 1 1];
            spec = sw_egrid(spinwave(sw_model('triAF', 1), {q0 q0+qdir testCase.nsteps{1}}));
            proj = [qdir(:) [1 -1 0]' [1 1 -2]'];
            dproj = [norm((qdir-q0))/testCase.nsteps{1}, 1e-6, 1e-6];
            sw_spec2MDHisto(spec, proj, dproj, 'tmp/test111mdh.nxs');
        end
        function test_1_1_2(testCase)
            q0 = [0 0 2];
            qdir = [1 1 0];
            spec = sw_egrid(spinwave(sw_model('triAF', 1), {q0 q0+qdir testCase.nsteps{1}}));
            proj = [qdir(:) [1 -1 0]' [0 0 1]'];
            dproj = [norm((qdir-q0))/testCase.nsteps{1}, 1e-6, 1e-6];
            sw_spec2MDHisto(spec, proj, dproj, 'tmp/test112mdh.nxs');
        end
        function test_1_1_2_2(testCase)
            q0 = [0 0 2];
            qdir = [1 1 0];
            spec = sw_egrid(spinwave(sw_model('triAF', 1), {q0 q0+qdir testCase.nsteps{1}}));
            proj = [[1 -1 0]' qdir(:)  [0 0 1]'];
            dproj = [1e-6, norm((qdir-q0))/testCase.nsteps{1}, 1e-6];
            sw_spec2MDHisto(spec, proj, dproj, 'tmp/test112_2mdh.nxs');
        end
        function test_2_2_2(testCase)
            q0 = [2 2 2];
            qdir = [1 1 0];
            spec = sw_egrid(spinwave(sw_model('triAF', 1), {q0 q0+qdir testCase.nsteps{1}}));
            proj = [qdir(:) [1 -1 0]' [0 0 1]'];
            dproj = [norm((qdir-q0))/testCase.nsteps{1}, 1e-6, 1e-6];
            sw_spec2MDHisto(spec, proj, dproj, 'tmp/test222mdh.nxs');
        end
        function test_non_ortho(testCase)
            q0 = [0 0 2];
            qdir = [1 1 0];
            spec = sw_egrid(spinwave(sw_model('triAF', 1), {q0 q0+qdir testCase.nsteps{1}}));
            proj = [qdir(:) [1 0 0]' [0 0 1]'];
            dproj = [1e-6, norm((qdir-q0))/testCase.nsteps{1}, 1e-6];
            verifyError(testCase,@() sw_spec2MDHisto(spec, proj, dproj, 'tmp/test_blank.nxs'), "read_struct:nonorthogonal")
        end
    end
end