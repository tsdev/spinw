%% setup help generator options

helpPath = {'swfiles/@spinw' 'swfiles' 'swfiles/+swplot' 'swfiles/+swpref' 'swfiles/+swsym' 'swfiles/+swfunc'};
swr      = sw_rootdir;
helpPath = cellfun(@(C)[swr C],helpPath,'UniformOutput',false);
swver    = sw_version;

%% generate help

%fun0 = {'spinw' 'spinw.spinwave' 'spinw.powspec'};
fun0 = cell(1,0);

clc

doctree = sw_genhelp('path',helpPath,'fun',fun0,'verstr',swver);


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



