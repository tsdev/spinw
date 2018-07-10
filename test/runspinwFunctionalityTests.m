function results = runspinwFunctionalityTests(d)
runner = matlab.unittest.TestRunner.withNoPlugins;
runner.addPlugin(matlab.unittest.plugins.TestSuiteProgressPlugin)
runner.addPlugin(matlab.unittest.plugins.FailureDiagnosticsPlugin)

if nargin == 0
    d = fullfile(sw_rootdir,'test','tutorialReport');
    if ~exist(d,'dir')
        mkdir(d);
    end
else
    if ~exist(d,'dir')
        error('spinw:FunctionalityTests', 'The specified directory %s is unavailable',d)
    end
end

tapFile = fullfile(d,'functionalityResults.tap');

runner.addPlugin(matlab.unittest.plugins.TAPPlugin.producingVersion13(matlab.unittest.plugins.ToFile(tapFile)))

results = runner.run(testsuite('spinwFunctionTests'));
end