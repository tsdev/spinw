function uiopen(type,direct)
% UIOPEN overloaded for custom Files. Do not change the file name of this
% file. Remember you are overloading uiopen inside toolbox/matlab/uitools
%
% You can use drag&drop for .cif files into the Matlab Command Window and
% they will be plotted automatically.
%

if ~isempty(strfind(type,'.cif')) && direct
    model = spinw(type);
    model.quickham([1 1 1]);
    plot(model);
else
    pwd0 = pwd;
    cd([matlabroot filesep 'toolbox' filesep 'matlab' filesep 'uitools'])
    feval('uiopen',type,direct);
    cd(pwd0)
end