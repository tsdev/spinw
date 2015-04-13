function gendoc_PSI(fName, DirName)
% generates documentation files for uploading onto the PSI SpinW page
%
% GENDOC_PSI({fName},{DirName})
%
%

if nargin < 2
    % list of files in the publish folder:
    pubfolder = [sw_rootdir 'tutorials' filesep 'help'];
else
    pubfolder = DirName;
end

if nargin == 0
    pubfiles  = dir([pubfolder filesep '*.m']);
else
    pubfiles  = dir([pubfolder filesep fName]);
end

% remove the publishing script
pubfiles(strcmp({pubfiles(:).name},'gendoc_PSI.m')) = [];

% open links.txt file and read expression list
fidL = fopen([pubfolder filesep 'links.txt']);
linklist = cell(0,2);
while ~feof(fidL)

    strL = fgets(fidL);
    linklist(end+1,:) = regexp(strL,'".*?"','match');
    linklist{end,1} = linklist{end,1}(2:end-1);
    linklist{end,2} = linklist{end,2}(2:end-1);
end
fclose(fidL);
%nLink = size(linklist,1);

% publish each file in separate folder
for ii = 1:numel(pubfiles)
    % exchange words with links
    fidIn  = fopen([pubfolder filesep pubfiles(ii).name]);
    fidOut = fopen([pubfolder filesep 'gen_' pubfiles(ii).name],'w');
    
    wordSel = find(~strcmpi(regexprep(linklist(:,1),'[^a-zA-Z]',''),pubfiles(ii).name(1:end-2)));
    
    while ~feof(fidIn)
        fLine = fgets(fidIn);
        %for jj = 1:nLink
        for jj = wordSel'
            % differentiate between sw.unit and sw.unit_cell
            fLine = regexprep(fLine,[linklist{jj,1} '(?=[^_A-Za-z])'],linklist{jj,2});
        end
        fprintf(fidOut,'%s',fLine);
    end
    
    fclose(fidOut);
    fclose(fidIn);
    
    use the non documented antialiasing option
    if verLessThan('Matlab','8.5')
        % doesn't work for R2015a or above
        opts.figureSnapMethod = 'antialiased';
    end
    opts.createThumbnail = false;
    opts.stylesheet = 'noCode.xsl';
    opts.outputDir = [pubfolder filesep pubfiles(ii).name(1:end-2)];
    publish([pubfolder filesep 'gen_' pubfiles(ii).name],opts);
    close all
    % keep only the text between <body> </body>
    fid = fopen([pubfolder filesep pubfiles(ii).name(1:end-2) filesep 'gen_' pubfiles(ii).name(1:end-1) 'html']);
    htmlFile = {};
    htmlFile{1} = fgets(fid);
    while ~feof(fid)
        htmlFile{end+1,1} = fgets(fid);
    end
    fclose(fid);
    
    % Matlab output default style format
    NL = sprintf('\n');
    htmlOut = {};
    htmlOut{end+1,1} = ['<literal>' NL];
    htmlOut{end+1,1} = ['<style type="text/css">' NL];
    htmlOut{end+1,1} = ['pre, tt, code { font-size:12px; }' NL];
    htmlOut{end+1,1} = ['pre.error { color:red; }' NL];
    htmlOut{end+1,1} = ['pre.codeinput  { padding:10px; border:1px solid #d3d3d3; background:#f7f7f7; }' NL];
    htmlOut{end+1,1} = ['pre.codeoutput { padding:10px; border:1px solid #d3d3d3; background:#FFFFFF; }' NL];
    htmlOut{end+1,1} = ['pre.language-matlab { padding:5px; border:0px solid #d3d3d3; background:#FFFFFF; }' NL];
    htmlOut{end+1,1} = ['span.keyword { color:#0000FF }' NL];
    htmlOut{end+1,1} = ['span.comment { color:#228B22 }' NL];
    htmlOut{end+1,1} = ['span.string { color:#A020F0 }' NL];
    htmlOut{end+1,1} = ['span.untermstring { color:#B20000 }' NL];
    htmlOut{end+1,1} = ['span.syscmd { color:#B28C00 }' NL];
    htmlOut{end+1,1} = ['</style>' NL];
    htmlOut{end+1,1} = ['</literal>' NL NL];

    % starting <body> index and line
    jj = 1;
    idx1 = strfind(htmlFile{jj},'<body>');
    while isempty(idx1)
        jj = jj+1;
        idx1 = strfind(htmlFile{jj},'<body>');
    end
    
    htmlFile{jj} = htmlFile{jj}(idx1+6:end);
    
    idx2 = [];
    while(isempty(idx2))
        
        htmlFile{jj} = regexprep(htmlFile{jj},'src="','src="%ATTACHURLPATH%/');
        idx2 = strfind(htmlFile{jj},'</body>');
        if ~isempty(idx2)
            htmlFile{jj} = htmlFile{jj}(1:idx2-1);
        end
        
        htmlOut{end+1,1} = htmlFile{jj}; %#ok<*AGROW>
        htmlOut{end,1} = regexprep(htmlOut{end,1},'LATEX','<latex>');
        htmlOut{end,1} = regexprep(htmlOut{end,1},'PATEX','</latex>');
        jj = jj+1;
    end
    
    % Google Analytics code
    htmlOut{end+1,1} = [NL NL '<literal>' NL];
    htmlOut{end+1,1} = ['<script>' NL];
    htmlOut{end+1,1} = ['(function(i,s,o,g,r,a,m){i[''GoogleAnalyticsObject'']=r;i[r]=i[r]||function(){' NL];
    htmlOut{end+1,1} = ['(i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),' NL];
    htmlOut{end+1,1} = ['m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)' NL];
    htmlOut{end+1,1} = ['})(window,document,''script'',''//www.google-analytics.com/analytics.js'',''ga'');' NL];
    htmlOut{end+1,1} = ['ga(''create'', ''UA-30000935-3'', ''psi.ch'');' NL];
    htmlOut{end+1,1} = ['ga(''send'', ''pageview'');' NL];
    htmlOut{end+1,1} = ['</script>' NL];
    htmlOut{end+1,1} = ['</literal>' NL];
    
    % create a single line
    htmlLine=sprintf('%s',htmlOut{:});
    
    fid = fopen([pubfolder filesep pubfiles(ii).name(1:end-2) filesep pubfiles(ii).name(1:end-1) 'txt'],'w');
    fprintf(fid,'%s',htmlLine);
    fclose(fid);
    delete([pubfolder filesep 'gen_' pubfiles(ii).name]);
    
end

end