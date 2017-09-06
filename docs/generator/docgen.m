%% generate help

swr = sw_rootdir;
clc
helpPath = {'swfiles' 'swfiles/@spinw' 'swfiles/+swplot' 'swfiles/+swpref' 'swfiles/+swsym'};
sw_genhelp('path',cellfun(@(C)[swr C],helpPath(4:5),'UniformOutput',false))

%%

