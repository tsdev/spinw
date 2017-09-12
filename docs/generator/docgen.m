%% generate help

swr = sw_rootdir;
clc
helpPath = {'swfiles' 'swfiles/@spinw' 'swfiles/+swplot' 'swfiles/+swpref' 'swfiles/+swsym'};
profile on
doctree = sw_genhelp('path',cellfun(@(C)[swr C],helpPath,'UniformOutput',false))
profile off

%%

