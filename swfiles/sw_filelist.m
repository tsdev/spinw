function varargout = sw_filelist(varargin)
% lists spinw data in the Matlab workspace or in a .mat file
%
% list = SW_FILELIST({'Option1',Value1,...)
%
%
% Options:
%
% fName     To check data stored in a .mat file, fName contains the path as
%           a string.
% sort      To sort according to which column (positive ascending, negative
%           descending):
%               +/-1    variable name,
%               +/-2    title,
%               +/-3    creation date,
%               +/-4    completion date, default.
%
%
% Output:
%
% list      Cell of strings, lists each simulation data in the Matlab
%           memory, even data stored in cells.
%
% See also SW.ANNEAL, SW.SPINWAVE.
%

inpForm.fname  = {'fName' 'sort' };
inpForm.defval = {''       +4    };
inpForm.size   = {[1 -1]  [1 1]  };
inpForm.soft   = {true    false  };

param = sw_readparam(inpForm, varargin{:});


if isempty(param.fName)
    S = evalin('base','whos');
    ws = 0;
else
    load(param.fName,'-mat');
    S = whos;
    ws = 1;
end

nVar = numel(S);

idx = 1;

varDateEnd    = {};
varDateStart  = {};
varTitle = {};
varName  = {};
idxDateEnd   = [];
idxDateStart = [];

oldDate = '01-Jan-1900 0:00:01';

for ii = 1:nVar
    if ws
        varS = eval(S(ii).name);
    else
        varS = evalin('base',S(ii).name);
    end
    
    if strcmp(S(ii).class,'struct') && isfield(varS,'obj')
        
        varName{idx} = S(ii).name; %#ok<*AGROW>
        % list variable
        if isfield(varS,'dateend')
            varDateEnd{idx}   = datestr(varS.dateend);
            varDateStart{idx} = datestr(varS.datestart);
            idxDateEnd(idx)   = datenum(varDateEnd{idx});
            idxDateStart(idx) = datenum(varDateStart{idx});
        else
            varDateEnd{idx}   = '-';
            varDateStart{idx} = '-';
            idxDateEnd(idx)   = datenum(oldDate);
            idxDateStart(idx) = datenum(oldDate);
        end
        
        if isfield(varS,'title')
            varTitle{idx} = varS.title;
        else
            varTitle{idx} = '-';
        end
        idx = idx + 1;
    end
    
    if strcmp(S(ii).class,'cell')
        if ws
            varC = eval(S(ii).name);
        else
            varC = evalin('base',S(ii).name);
        end
        for jj = 1:numel(varC)
            if isfield(varC{jj},'obj')
                % list variable from cell
                
                varName{idx} = [S(ii).name '{' num2str(jj) '}'];
                if isfield(varC{jj},'dateend')
                    varDateEnd{idx}   = datestr(varC{jj}.dateend);
                    varDateStart{idx} = datestr(varC{jj}.datestart);
                    idxDateEnd(idx)   = datenum(varDateEnd{idx});
                    idxDateStart(idx) = datenum(varDateStart{idx});
                else
                    varDateEnd{idx}   = '-';
                    varDateStart{idx} = '-';
                    idxDateEnd(idx)   = datenum(oldDate);
                    idxDateStart(idx) = datenum(oldDate);
                end
                
                if isfield(varC{jj},'title')
                    varTitle{idx} = varC{jj}.title;
                else
                    varTitle{idx} = '-';
                end
                idx = idx + 1;
                
            end
        end
    end
end

% sort data according to selection
switch param.sort
    case {-1, 1}
        [~, idxs] = sort(varName);
    case {-2, 2}
        [~, idxs] = sort(varTitle);
    case {-3, 3}
        [~, idxs] = sort(idxDateStart);
    case {-4, 4}
        [~, idxs] = sort(idxDateEnd);
    otherwise
end

if sign(param.sort)<0
    idxs = idxs(end:-1:1);
end

idxDateEnd = idxDateEnd(idxs);
varDateEnd   = varDateEnd(idxs);
varDateStart = varDateStart(idxs);
varName    = varName(idxs);
varTitle   = varTitle(idxs);

% no variable found
if idx < 2
    if isempty(param.fName)
        fprintf('No SpinW variables have been found in the Matlab base workspace.\n');
    else
        fprintf('No SpinW variables stored in the "%s" file.\n',param.fName);
    end
    return
end

if nargout == 0
    if isempty(param.fName)
        fprintf('SpinW variables in the Matlab base workspace:\n');
    else
        fprintf('SpinW variables stored in the "%s" file:\n',param.fName);
    end
    fprintf('%15s %45s %25s %25s\n','VarName','Title','Creation Date','Completion Date');
    for ii = 1:idx-1
        strSp = strsplit(varName{ii},'{');
        if numel(strSp)==1
            strSp{2} = '';
        else
            strSp{2} = ['{' strSp{2}];
        end
        varTitle0 = varTitle{ii};
        if numel(varTitle0)>37
            varTitle0 = [varTitle0(1:37) '...'];
        end
        fprintf('%15s%-5s %40s %25s %25s\n',strSp{:},varTitle0,varDateStart{ii},varDateEnd{ii})
    end
else
    varargout{1}.name = varName;
    varargout{1}.dateend   = varDateEnd;
    varargout{1}.datestart = varDateStart;
    varargout{1}.datenum = idxDateEnd;
    varargout{1}.title = varTitle;
end

end