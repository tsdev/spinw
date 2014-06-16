function publish_PSI(fName)
% publishing all examples in the publish folder and convert it to PSI web format
%
% The online source can be imported to Matlab using the web address of the
% tutorial file:
%
% Example:
% grabcode('http://www.psi.ch/spinw/tutorial-2')
%
% creates .txt files for publishing on the PSI website
%

% list of files in the publish folder:
pubfolder = [sw_rootdir 'doc' filesep 'publish'];
if nargin == 0
    pubfiles  = dir([pubfolder filesep '*.m']);
else
    pubfiles  = dir([pubfolder filesep fName]);
end

% remove the publishing script
pubfiles(strcmp({pubfiles(:).name},'publish_PSI.m')) = [];


% publish each file in separate folder
for ii = 1:numel(pubfiles)
    %publish([pubfolder filesep pubfiles(ii).name],'outputDir',[pubfolder filesep pubfiles(ii).name(1:end-2)],'maxOutputLines',0);
    % use the non documented antialiasing option
    opts.figureSnapMethod = 'antialiased';
    opts.outputDir = [pubfolder filesep pubfiles(ii).name(1:end-2)];
    publish([pubfolder filesep pubfiles(ii).name],opts);
    close all
    % keep only the text between <body> </body>
    fid = fopen([pubfolder filesep pubfiles(ii).name(1:end-2) filesep pubfiles(ii).name(1:end-1) 'html']);
    htmlFile = {};
    htmlFile{1} = fgets(fid);
    while ~feof(fid)
        htmlFile{end+1,1} = fgets(fid);
    end
    fclose(fid);
    
    htmlOut = {};
    % starting <body> index and line
    jj = 1;
    idx1 = strfind(htmlFile{jj},'<body>');
    while isempty(idx1)
        jj = jj+1;
        idx1 = strfind(htmlFile{jj},'<body>');
    end
    htmlOut{1} = htmlFile{jj}(idx1+6:end);
    jj = jj + 1;
    
    % ending </body> index and line
    idx2 = strfind(htmlFile{jj},'</body>');
    imgPath = '%ATTACHURLPATH%/';
    while(isempty(idx2))
        % change image path
        imIdx = strfind(htmlFile{jj},'src="');
        if ~isempty(imIdx)
            lineSel = htmlFile{jj};
            
            %             switch numel(imIdx)
            %                 case 1
            %                     htmlFile{jj} = [lineSel(1:imIdx+4) imgPath lineSel(imIdx+5:end)];
            %                 case 2
            %                     htmlFile{jj} = [lineSel(1:imIdx(1)+4) imgPath lineSel(imIdx(1)+5:imIdx(2)+4) imgPath lineSel(imIdx(2)+5:end)];
            %                 case 3
            %                     htmlFile{jj} = [lineSel(1:imIdx(1)+4) imgPath lineSel(imIdx(1)+5:imIdx(2)+4) imgPath lineSel(imIdx(2)+5:imIdx(3)+4) imgPath lineSel(imIdx(3)+5:end)];
            %                 otherwise
            %             end
            
            htmlFile{jj} = [lineSel(1:imIdx(1)+4) imgPath];
            idxSum = 1;
            while idxSum < numel(imIdx)
                htmlFile{jj} = [htmlFile{jj} lineSel(imIdx(idxSum)+5:imIdx(idxSum+1)+4) imgPath];
                idxSum = idxSum + 1;
            end
            htmlFile{jj} = [htmlFile{jj} lineSel(imIdx(end)+5:end)];
            
        end
        htmlOut{end+1,1} = htmlFile{jj}; %#ok<*AGROW>
        jj = jj+1;
        idx2 = strfind(htmlFile{jj},'</body>');
    end
    htmlOut{end+1,1} = htmlFile{jj}(1:idx2-1);
    % create a single line
    htmlLine=sprintf('%s',htmlOut{:});
    %htmlLine = [];
    %for jj = 1:numel(htmlOut)
    %    htmlLine = [htmlLine htmlOut{jj} '\n'];
    %end
    
    % put comment for the PSI publishing script to keep original html code
    srcStart = regexp(htmlLine,'##### SOURCE BEGIN #####');
    srcEnd   = regexp(htmlLine,'##### SOURCE END #####');
    htmlLine = [htmlLine(1:srcStart-1) '<literal>' htmlLine(srcStart:srcEnd+21) '</literal>' htmlLine(srcEnd+22:end)];
    
    fid = fopen([pubfolder filesep pubfiles(ii).name(1:end-2) filesep pubfiles(ii).name(1:end-1) 'txt'],'w');
    %for jj = 1:numel(htmlOut)
    fprintf(fid,'%s',htmlLine);
    %end
    fclose(fid);
    
end

end