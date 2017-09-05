function sw_genhelp(varargin)
% generates markdown files from function help
%
% SW_GENHELP('option1', value1, ...)
%

inpForm.fname  = {'path'};
inpForm.defval = {{}    };
inpForm.size   = {[1 -1]};

param = sw_readparam(inpForm, varargin{:});

if ~iscell(param.path)
    path0 = {param.path};
else
    path0 = param.path;
end

docroot = [sw_rootdir 'docs' filesep];

% loop over all path to generate help files
for ii = 1:numel(path0)
    % name of the parent folder
    [~,pp1,pp2] = fileparts(path0{ii});
    folder = [pp1 pp2];
    
    % find all *.m files in the folder
    fList = dir([path0{ii} filesep '*.m']);
    fList = {fList(:).name};
    
    % remove old help
    try
        rmdir([docroot 'pages' filesep folder],'s')
    catch
    end
    % add new folder
    mkdir([docroot 'pages' filesep folder]);
    % load the help file for each file and save them as a markdown file
    
    for jj = 1:numel(fList)
        helpText = help([path0{ii} filesep fList{jj}]);
        % generate title
        summary = strsplit(helpText,'\n');
        summary = strtrim(summary{1});
        funName = fList{jj}(1:end-2);
        header = ['---\ntitle: ' funName '( )\nkeywords: sample\nsummary: "' summary '"\nsidebar: product1_sidebar\npermalink: ' funName '.html\nfolder: ' folder '\nmathjax: true\n---\n'];
        
        % save the text as .md file
        fid = fopen([docroot 'pages' filesep folder filesep funName '.md'],'w');
        fprintf(fid,[header helpText]);
        fclose(fid);
    end
    
end
















end