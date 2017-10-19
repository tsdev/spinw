%% setup help generator options

swPath  = {'swfiles/@spinw' 'swfiles' 'swfiles/+swplot' 'swfiles/+swpref' 'swfiles/+swsym' 'swfiles/+swfunc'};
swr     = sw_rootdir;
swPath  = cellfun(@(C)[swr C],swPath,'UniformOutput',false);
swver   = sw_version;
outPath = '~/spinwdoc_git';
docPath = '~/spinw_git/docs';
upload  = false;
recalc  = false;

%% generate help

%fun0 = {'swfiles' 'spinw' 'gm_planard' 'sw_version'};
%fun0 = {'swfiles' 'sw_egrid'};
fun0 = cell(1,0);

clc

doctree = sw_genhelp('sourcepath',swPath,'outpath',outPath,'docpath',docPath,'fun',fun0,'verstr',swver,'recalc',recalc);

% upload git
if upload
    pwd0 = pwd;
    cd(outPath)
    !git add .
    !git commit -m "Regenerated documentation"
    !git push
    cd(pwd0)
end
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



