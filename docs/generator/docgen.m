%% generate help

swr = sw_rootdir;
clc
helpPath = {'swfiles' 'swfiles/@spinw' 'swfiles/+swplot' 'swfiles/+swpref' 'swfiles/+swsym' 'swfiles/+swfunc'};
profile on
doctree = sw_genhelp('path',cellfun(@(C)[swr C],helpPath,'UniformOutput',false));
profile off

%% get all help

content = [doctree.content];
allhelp = {content.text};

for ii = 1:numel(allhelp)
    allhelp{ii} = cellfun(@(C)[C newline],allhelp{ii},'UniformOutput',false);
    allhelp{ii} = [allhelp{ii}{:}];
end

%%

tic
for ii = 1:215
    
    regexprep(allhelp,'version[ ]+-','[spinw.version] -');
end
toc