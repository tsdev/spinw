classdef spinwFunctionTests < matlab.unittest.TestCase
    
    
    properties
        generate_ref = false;
        seed = 3142;
        relToll = 0.025;
        absToll = 1E-6;
        refDir = fullfile(sw_rootdir, 'test', 'tutorialResuls');
    end
    
    properties(Access = protected)
       r1 
    end
    
    properties (TestParameter)
        testFun = arrayfun(@(x) fullfile(x.folder, x.name), ...
            dir(fullfile(sw_rootdir, 'tutorials', 'publish', 'tutorial*.m')),...
            'UniformOutput', false);
%         Incase we need to debug one tutorial.
%         testFun = {fullfile(sw_rootdir, 'tutorials', 'publish', 'tutorial9.m')};
        
    end
    
    methods (TestMethodSetup)
        function prepareForRun(testCase)
            
           close all
           clc
           
           testCase.r1 = rng;
           rng(testCase.seed)
           
        end
    end
    
    methods (TestMethodTeardown)
        function clearnupRun(testCase)
           rng(testCase.r1)
        end
    end
    
    methods (TestClassSetup)
        function checkDirs(testCase)
            if ~exist(testCase.refDir,'dir')
               testCase.generate_ref = true;
               warning('Making references as none exist')
            end
        end
    end
    
    methods (TestClassTeardown)
        function closeAllFigs(testCase)
           close all
        end
    end
    
    methods (Test)
        function testFunction(testCase, testFun)
            outStat = 1;
            
            str= fileread(which(testFun));  %read in text file to memory
            urls = regexp(str, 'http(\S+)(\s*)$', 'match', 'lineAnchors');    %find urls
            
            if ~isempty(urls)
                for i = 1:length(urls)
                    url = strtok(urls{i},''')');
                    if ~contains(url,'goo.gl')
                        break
                    end
                    try
                        [~] = webread(strtok(url,''')'));
                    catch ME
                        ident = '';
                        switch ME.identifier
                            case 'MATLAB:webservices:HTTP404StatusCodeError'
                                warning('WARNING! Remote content not found for tutorial %s', testFun);
                                ident = 'MATLAB:webservices:HTTP404StatusCodeError';
                            case 'MATLAB:webservices:UnknownHost'
                                warning('WARNING! Can not find remote content. DNS problem?')
                                ident = 'MATLAB:webservices:UnknownHost';
                            case 'MATLAB:webservices:Timeout'
                                warning('WARNING! Can not connect to the internet.')
                                ident = 'MATLAB:webservices:Timeout';
                            otherwise
                                rethrow(ME)
                        end
                        import matlab.unittest.constraints.Matches
                        testCase.verifyMatches(ME.identifier, ident);
                        return
                    end
                end
            end
            
            try
                tic
                run(testFun)
                t = toc;
            catch ME
                outStat = 0;
                sprintf('%s',ME.message)
            end
            
            % Fail if the test fails
            testCase.verifyEqual(outStat, 1);
            
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            
            theseBounds = RelativeTolerance(testCase.relToll) | AbsoluteTolerance(testCase.absToll);
            
            % Verify Generated data is accurate
            axes = findobj('type','axes');
            axes = axes(arrayfun(@(x) isempty(x.Tag),axes));
                        
            for i = 1:length(axes)
                data_s = extractData(axes(i));
                temp = strtok(testFun(end:-1:1), filesep);
                r_name = fullfile(testCase.refDir, sprintf('%s_ref%i.mat', temp(end:-1:3), i));
                if testCase.generate_ref
                    data_r = data_s;
                    data_r = arrayfun(@(x) setfield(x, 'matVer', version), data_r);
                    data_r = arrayfun(@(x) setfield(x, 'time', t), data_r);
                    data_r = arrayfun(@(x) setfield(x, 'spinwVer', sw_version), data_r);
                    save(r_name, 'data_r')
                else
                    load(r_name, 'data_r');
                    for j = 1:length(data_s)
                        testCase.verifyThat(data_s(j).data, IsEqualTo(data_r(j).data, 'Within', theseBounds))
                    end
                    clear data_r;
                end
            end
        end
    end
end

function dataOut = extractData(axisIN)
children = axisIN.Children;
for i = 1:length(children)
    child = children(i);
    dataOut(i).axisNo = i;
    switch get(child,'Type')
        case 'line'
            dataOut(i).data = [child.XData(:) child.YData(:)];
        case 'surface'
            dataOut(i).data = [child.XData(:) child.YData(:) child.ZData(:)];
        otherwise
            dataOut(i).data = [];
    end
end
end