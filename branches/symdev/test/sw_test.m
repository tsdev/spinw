function errMsg = sw_test(varargin)
% testing of sw functions
% errMsg = sw_test('Option1,' Value1, ...)
%
% Options:
%
% fid       Identifier where the text output goes.
%           0   No output.
%           1   Output onto the Command Window.
%           fid Output into file, opened with: fid = fopen(path)
% tol       Tolerance on the agreement of different calculated matrices,
%           default is 1e-5.
%
% Output:
%
% errMsg    Cell, that contains all the error messages for every test
%           function.
%

currDir = cd;
cd(sw_rootdir);

% switch off warnings
warnLevel = warning;
warning('off','all');

inpForm.fname  = {'fid' 'tol' };
inpForm.defval = {1     1e-5  };
inpForm.size   = {[1 1] [1 1] };

param = sw_readparam(inpForm, varargin{:});

fid = param.fid;

Res = [];
errMsg = {};

allTestPath = dir([sw_rootdir 'test' filesep 'sw_test_*.m']);
allTestPath = {allTestPath.name};

for ii = 1:numel(allTestPath)
    [Res(ii), errMsg{ii}] = eval([allTestPath{ii}(1:end-2) '(param.tol)']); %#ok<AGROW>
end

minRes = min(Res(Res>0));
if isempty(minRes)
    minRes = 0;
end

switch minRes
    case 0
        % no error
        fprintf(fid,'All test ran succesfull!\n');
    case 1
        % error message
        fprintf(fid,'The code throw an error message!\n');
    case 2
        % wrong result
        fprintf(fid,'Some of the numerical results were wrong!\n');
    otherwise
        fprintf(fid,'Unkown error code!\n');
end

errMsg = [allTestPath; errMsg];

% clear symmetry.dat file from newly added entries
sw_initialize;

% change back to the current directory
cd(currDir);

% switch on warnings
warning(warnLevel);

end