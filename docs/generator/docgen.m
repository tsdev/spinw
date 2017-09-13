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
isCont  = [content.isContents];

for ii = 1:numel(allhelp)
    allhelp{ii} = cellfun(@(C)[C newline],allhelp{ii},'UniformOutput',false);
    allhelp{ii} = [allhelp{ii}{:}];
end

%% generate all links

pLink = cellfun(@(C)C.permalink,{content.frontmatter},'UniformOutput',false);
fun   = {content.fun};

%% replace function names with links

for ii = 1:numel(pLink)
    allhelp = regexprep(allhelp,fun{ii},['[' fun{ii} '](' pLink{ii} ')']);
end



