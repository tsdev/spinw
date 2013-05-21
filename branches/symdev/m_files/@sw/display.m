function display(obj, fid)
%  prints the sw data structure in readable format onto the Command Window
%

if nargin == 1
    fid = 1;
end

Datastruct = datastruct;
fprintf(fid, 'sw object\n');

for ii = 1:length(Datastruct.mainfield)
    fprintf(fid, '%s\n', Datastruct.mainfield{ii});
    index = 1;
    plotSize = true;
    while index <= size(Datastruct.subfield,2) && ~isempty(Datastruct.subfield{ii,index})
        fprintf(fid, '    %10s: ',Datastruct.subfield{ii,index});
        sizefield = Datastruct.sizefield{ii,index};
        typefield = Datastruct.typefield{ii,index};
        
        fprintf(fid, '[');
        
        for jj = 1:length(sizefield)-1
            objSize = [size(obj.(Datastruct.mainfield{ii}).(Datastruct.subfield{ii,index})) 1 1 1];
            if ischar(sizefield{jj})
                if plotSize
                    %fprintf(fid, '%s=%dx', sizefield{jj},objSize(jj));
                    fprintf(fid, '%sx', sizefield{jj});
                    sF = sizefield{jj};
                    sS = objSize(jj);
                else
                    fprintf(fid, '%sx', sizefield{jj});
                end
                
            else
                fprintf(fid, '%dx', sizefield{jj});
            end
        end
        
        if ischar(sizefield{jj+1})
            if plotSize
                fprintf(fid, '%s', sizefield{jj+1});
                sF = sizefield{jj+1};
                sS = objSize(jj+1);
                
            else
                fprintf(fid, '%s', sizefield{jj+1});
            end
            
        else
            fprintf(fid, '%d', sizefield{jj+1});
        end
        fprintf(fid, ' %s]',typefield);
        if plotSize && exist('sF','var')
            fprintf('  %s=%d',sF,sS);
            clear sF sS
            plotSize = false;
        end
        fprintf('\n')
        index = index + 1;
    end
end

end