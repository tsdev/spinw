function cifdat = importcif(path)
% imports .cif file

fid = fopen(path);

loop = false;
loopvar = {};


idx = 1;
cifStr = {};

% read in all lines
while ~feof(fid)
    cifStr{idx} = fgetl(fid); %#ok<AGROW>
    idx = idx + 1;
end
% unite broken lines
bLine = cellfun(@(x)x(1)==';',cifStr);

% first and last line of the series of broken lines
firstBL = find(diff([0 bLine])>0);
lastBL  = find(diff([bLine 0])<0);

cifStr2 = {};
BL = false;
idx = 1;

for ii = 1:numel(cifStr)
    if any(ii==firstBL)
        BL = true;
    end
    
    if ~BL
        cifStr2{end+1} = cifStr{ii};
    else
        % end comment sign
        cifStr2{end} = [cifStr2{end} ' ±'  cifStr{ii}(2:end)];
    end
    
    if any(ii==lastBL)
        BL = false;
    end
end

cifStr = cifStr2;

strout = {};
for ii = 1:numel(cifStr)
    strout{ii} = parseline(cifStr{ii});
end

% remove comments
for ii = 1:numel(strout)
    comIdx = strcmp('comment',strout{ii}(1,:));
    strout{ii} = strout{ii}(:,~comIdx);
end

% remove empty lines
emptyIdx = cellfun(@(x)isempty(x),strout);
strout(emptyIdx) = [];

LP = false;
lpVar = {};
lpVal = {};

cifdat.name = '';
cifdat.val  = '';

for ii = 1:numel(strout)
    strLine = strout{ii};
    switch strLine{1,1}
        case 'loop'
            LP = true;
            lpVar = {};
            lpVal = {};
        case 'variable'
            if LP
                lpVar = [lpVar strLine(2,1)];
            elseif (size(strLine,2)> 1) && ismember(strLine(1,2),{'number' 'string'})
                cifdat(end+1).name = strLine{2,1};
                cifdat(end).val    = strLine{2,2};
            else
                cifdat(end+1).name = strLine{2,1};
                cifdat(end).val    = '';
                warning('Wrong stuff comes after a variable in line: %d',ii);
            end
            
            % end loop
            if ~isempty(lpVal)
                LP = false;
                % save data
                for jj = 1:numel(lpVar)
                    cifdat(end+1).name = lpVar{jj};
                    cifdat(end).val    = lpVal(:,jj);
                end
                lpVar = {};
                lpVal = {};
                
            end
        case 'number'
            if LP
                lpVal = {lpVal; strLine(2,:)};
            end
        case 'string'
            if LP
                lpVal = {lpVal; strLine(2,:)};
            end
            
        otherwise
            if ~isempty(lpVal)
                % save data
                for jj = 1:numel(lpVar)
                    cifdat(end+1).name = lpVar{jj};
                    cifdat(end).val    = lpVal(:,jj);
                end
                lpVar = {};
                lpVal = {};
            end
            LP = false;
    end
end


fclose(fid);
cifdat = cifdat(2:end);

end



function strout = parseline(strin)
% parse cif file line
strout = cell(2,0);

while ~isempty(strin)
    % remove whitespace
    strin = strtrim(strin);
    
    if isempty(strin)
        return;
    end
    
    % COMMENT
    if strin(1) == '#'
        endcomIdx    = find(strin(1:end)=='±',1,'first');
        if isempty(endcomIdx)
            endcomIdx = numel(strin)+1;
        end
        
        strout{1,end+1} = 'comment';
        strout{2,end}   = strin(2:(endcomIdx-1));
        strin = strin(endcomIdx:end);
        continue;
    end
    
    % VARIABLE
    if strin(1) == '_'
        whiteIdx    = find(strin(1:end)==' ',1,'first');
        if isempty(whiteIdx)
            whiteIdx = numel(strin)+1;
        end
        strout{1,end+1} = 'variable';
        strout{2,end} = strin(2:(whiteIdx-1));
        strin = strin((whiteIdx+1):end);
        continue;
    end
    
    % LOOP
    if (numel(strin) >4) && strcmpi(strin(1:5),'loop_')
        strout{1,end+1} = 'loop';
        strout{2,end} = [];
        return;
    end
    
    % STRING
    if strin(1) == ''''
        strIdx = find(strin(2:end)=='''',1,'first')+1;
        if isempty(strIdx)
            error('importcif:WrongCif','String is not closed!');
        end
        
        strout{1,end+1} = 'string';
        strout{2,end} = strin(2:(strIdx-1));
        strin = strin((strIdx+1):end);
        continue;
    end
    
    % NUMBER
    if ismember(strin(1), [mat2cell('0':'9',1,ones(1,10)) {'-'}])
        bracket1Idx = find(strin(1:end)=='(',1,'first');
        bracket2Idx = find(strin(1:end)==')',1,'first');
        whiteIdx    = find(strin(1:end)==' ',1,'first');
        if isempty(bracket1Idx)
            bracket1Idx = numel(strin)+1;
        end
        if isempty(bracket2Idx)
            bracket2Idx = numel(strin)+1;
        end
        if isempty(whiteIdx)
            whiteIdx = numel(strin)+1;
        end
        
        if whiteIdx <= bracket1Idx
            strout{1,end+1} = 'number';
            strout{2,end} = str2num(strin(1:(whiteIdx-1)));
            strin = strin((whiteIdx+1):end);
            continue;
        else
            strout{1,end+1} = 'number';
            strout{2,end} = str2num(strin(1:(bracket1Idx-1)));
            strin = strin((bracket2Idx+1):end);
            continue;
            
        end
    end
    
    % STRING
    if isstrprop(strin(1), 'alpha')
        whiteIdx    = find(strin(1:end)==' ',1,'first');
        if isempty(whiteIdx)
            whiteIdx = numel(strin)+1;
        end
        
        strout{1,end+1} = 'string';
        strout{2,end} = strin(1:(whiteIdx-1));
        strin = strin((whiteIdx+1):end);
        continue;
    end
    
    % EMPTY NUMBER
    if strin(1) == '?'
        whiteIdx    = find(strin(1:end)==' ',1,'first');
        strout{1,end+1} = 'number';
        strout{2,end} = [];
        strin = strin((whiteIdx+1):end);
        continue;
    end
    
    % END COMMENT WITHOUT BEGIN
    if strin(1) == '±'
        strin = strin(2:end);
        continue;
    end
    
end


end

