function [str, out] = sw_example(str,path,figname,recalc)
% executes example code given in a cell of strings
%
% Lines starting with >> will be executed
% Lines starting with >>> will be executed, but not shown
% Line containing snapnow causes to save the image of the last figure

keep = true(1,numel(str));

out = {};

warn0 = warning;
warning off
cc = onCleanup(@() warning(warn0));

% find code blocks
idx = 1;

for ii = 1:numel(str)
    temp = strtrim(str{ii});
    if ~isempty(regexp(temp,'snapnow','once'))
        % save figure
        fignamei = [figname '_' num2str(idx) '.png'];
        
        if recalc    
            idx = idx+1;
            set(gcf,'Color',[241 240 240]/255*0+1)
            drawnow;
            ihc = get(gcf,'InvertHardCopy');
            set(gcf,'InvertHardCopy','off');
            print([path filesep fignamei],'-dpng','-r144');
            set(gcf,'InvertHardCopy',ihc);
            close(gcf);
        end
        caption = str{ii-1};
        caption = strtrim(caption(caption~='>'));
        % removed --> max-width="500"
        str{ii} = ['```' newline ' ' newline '{% include image.html file="' fignamei '" alt="' caption '" %}' newline newline '```matlab' newline];
        keep(ii) = true;
    elseif numel(temp)>3 && strcmp(temp(1:3),'>>>')
        if recalc
            evalc(str{ii}(4:end));
        end
        keep(ii) = false;
    elseif numel(temp)>2 && strcmp(temp(1:2),'>>')
        if recalc
            out0 = evalc(str{ii}(str{ii}~='>'));
            if strcmp(temp(end+(-1:0)),'>>')
                % add output to the string, remove html tags
                out1 = regexprep(out0,'<.*?>','');
                out1 = strsplit(out1,newline);
                out1 = out1(cellfun(@(C)isempty(C),regexp(out1,'^ans')));
                out1 = sprintf('%s\n',out1{:});
                str{ii} = [str{ii} newline '```' newline '*Output*' newline '```' newline out1(1:end-1) newline '```' newline ' ' newline '```matlab'];
            else
                % just keep it in out
                
            end
        else
            out0 = '';
        end
        if ~isempty(out0)
            out{end+1} = out0; %#ok<AGROW>
        end
        str{ii} = str{ii}(str{ii}~='>');
    end
end

% close all opened figures
close('all')
str = str(keep);
str = strsplit(sprintf('%s\n',str{:}),newline)';
str = str(1:end-1);


end