function varargout = display(obj, fid)
% prints the sw data structure in readable format onto the Command Window
%
% {swDescr} = DISPLAY(obj, {fid})
%
% Input:
%
% obj       sw class object.
% fid       File identifier (created by fopen function) where the output
%           will be written. Default is 1, the output will be written onto
%           the Command Line. Optional
%
% Output:
%
% swdescr   If output variable is given, the description of the obj object
%           will be output into the swdescr variable, instead of being
%           written onto the Command Window/file. Optional.
%
% Example:
%
% crystal = sw;
% swFields = display(crystal);
%
% See also SW.
%

Datastruct = datastruct;
swDescr = sprintf('sw object\n');

for ii = 1:length(Datastruct.mainfield)
    swDescr = [swDescr sprintf('%s\n', Datastruct.mainfield{ii})]; %#ok<*AGROW>
    index = 1;
    plotSize = true;
    while index <= size(Datastruct.subfield,2) && ~isempty(Datastruct.subfield{ii,index})
        swDescr = [swDescr sprintf('    %10s: ',Datastruct.subfield{ii,index})];
        sizefield = Datastruct.sizefield{ii,index};
        typefield = Datastruct.typefield{ii,index};
        
        swDescr = [swDescr sprintf('[')];
        
        for jj = 1:length(sizefield)-1
            objSize = [size(obj.(Datastruct.mainfield{ii}).(Datastruct.subfield{ii,index})) 1 1 1];
            if ischar(sizefield{jj})
                if plotSize
                    swDescr = [swDescr sprintf('%sx', sizefield{jj})];
                    sF = sizefield{jj};
                    sS = objSize(jj);
                else
                    swDescr = [swDescr sprintf('%sx', sizefield{jj})];
                end
                
            else
                swDescr = [swDescr sprintf('%dx', sizefield{jj})];
            end
        end
        
        if ischar(sizefield{jj+1})
            if plotSize
                swDescr = [swDescr sprintf('%s', sizefield{jj+1})];
                sF = sizefield{jj+1};
                sS = objSize(jj+1);
                
            else
                swDescr = [swDescr sprintf('%s', sizefield{jj+1})];
            end
            
        else
            swDescr = [swDescr sprintf('%d', sizefield{jj+1})];
        end
        swDescr = [swDescr sprintf(' %s]',typefield)];
        if plotSize && exist('sF','var')
            swDescr = [swDescr sprintf('  %s=%d',sF,sS)];
            clear sF sS
            plotSize = false;
        end
        swDescr = [swDescr sprintf('\n')];
        index = index + 1;
    end
end

if nargout == 1
    varargout{1} = swDescr;
else
    if nargin == 1
        fid = 1;
    end
    % print the text
    fprintf(fid,swDescr);
end

end