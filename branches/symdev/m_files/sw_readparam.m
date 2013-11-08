function input = sw_readparam(format, varargin)
% input = SW_READPARAM(format, raw) reads in parameters from input
% structure. Lower and upper case insensitive, the output structure has the
% names stored in format.fname. Instead of a struct type input, also a list
% of parmeters can be given in a parameter name, value pairs. Where the
% parameter name is a string.
%
% Input:
%
% format is struct type with the following fields:
% fname     Field names, strings in cell, dimensions are [nParm 1].
% size      field size, if negative means index, field sizes with same
%           negative index have to be the same size.
% defval    Optional, default value if missing.
% soft      Optional, if exist and equal to 1, in case of bad input
%           value, defval is used without error message.
%

if nargin == 0
    help sw_readparam;
    return;
end

% create a showWarn field to check whether to show warnings (default true)
format.fname  = [format.fname  {'showWarn'}];
format.defval = [format.defval {true      }];
format.size   = [format.size   {[1 1]     }];
if isfield(format,'soft')
    format.soft   = [format.soft   {true      }];
end

if (nargin>2) && (mod(nargin,2) == 1)
    nPar = nargin-1;
    raw = struct;
    for ii = 1:2:nPar
        raw.(varargin{ii}) = varargin{ii+1};
    end
elseif nargin == 2
    raw = varargin{1};
elseif nargin == 1
    raw = struct;
else
    error('sw:sw_readparam:WrongNumberOfInput','Wrong number of input parameters!');
end

fName     = format.fname;
rName     = fieldnames(raw);
storeSize = zeros(20,1);
input     = struct;

usedField = false(1,numel(rName));

% Go through all fields.
for ii = 1:length(fName)
    
    rawIdx = find(strcmpi(rName,fName{ii}));
    
    if any(rawIdx)
        rawIdx = rawIdx(1);
        usedField(rawIdx) = true;
        
        inputValid = true;
        
        % Go through all dimension of the selected field to check size.
        for jj = 1:length(format.size{ii})
            if format.size{ii}(jj)>0
                if format.size{ii}(jj) ~= size(raw.(rName{rawIdx}),jj)
                    inputValid = false;
                end
            else
                if storeSize(-format.size{ii}(jj)) == 0
                    storeSize(-format.size{ii}(jj)) = size(raw.(rName{rawIdx}),jj);
                else
                    if storeSize(-format.size{ii}(jj)) ~= size(raw.(rName{rawIdx}),jj)
                        inputValid = false;
                    end
                    
                end
            end
        end
        
        if inputValid
            input.(fName{ii}) = raw.(rName{rawIdx});
        else
            if isfield(format,'soft') && format.soft{ii}
                input.(fName{ii}) = format.defval{ii};
            else
                error('spinw:sw_readparam:SizeMismatch',['Input parameter size mismatch in .' fName{ii} '!']);
            end
        end
    else
        if isfield(format,'defval') && (any(size(format.defval{ii})) || (isfield(format,'soft') && format.soft{ii}))
            input.(fName{ii}) = format.defval{ii};
        else
            error('spinw:sw_readparam:ParameterMissing',['Necessary input parameter missing (param.' fName{ii} ')!']);
        end
    end
end

if input.showWarn && ~all(usedField)
    wName = sprintf('%s, ',rName{~usedField});
    warning('sw_readparam:UnreadInput','Invalid input parameter names: %s\b\b!',wName);
end

input = rmfield(input,'showWarn');

end