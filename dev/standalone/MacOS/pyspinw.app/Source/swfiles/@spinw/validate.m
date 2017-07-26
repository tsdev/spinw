function validate(varargin)
% validates spinw object properties
%
% VALIDATE(obj, {fieldToValidate})
%

% Load data structure.
Datstruct = datastruct();

mainfield = Datstruct.mainfield;
subfield  = Datstruct.subfield;
sizefield = Datstruct.sizefield;
typefield = Datstruct.typefield;



% Validate only selected mainfield of the struct.
indexFieldToValidate = [];
if nargin > 1
    if isa(varargin{2},'char')
        fieldToValidate = varargin(2);
    end
    
    for ii = 1:length(fieldToValidate)
        indexFieldToValidate = [indexFieldToValidate find(strcmp(mainfield,fieldToValidate{ii}))]; %#ok<AGROW>
        
    end
end

if isempty(indexFieldToValidate)
    indexFieldToValidate = 1:length(mainfield);
end

valid  = true;
objS   = struct(varargin{1});
fieldM = '';

for ii = indexFieldToValidate
    selectMainField = mainfield{ii};
    validT = isfield(objS,selectMainField);
    if valid && ~validT
        fieldM = selectMainField;
    end
    valid = valid && validT;
    
    for jj = 1:size(subfield,2)
        selectSubField = subfield{ii,jj};
        if ~isempty(selectSubField)
            validT = isfield(objS.(selectMainField),selectSubField);
            if valid && ~validT
                fieldM = [selectMainField '.' selectSubField];
            end
            valid  = valid && validT;
            
        end
    end
end

if ~valid
    error('spinw:sw_valid:MissingField',['Input struct missing necessary field: ' fieldM '!'])
end

% Save the type of error for easier debugging.
errorType = 0;
errorData = '';

for ii = indexFieldToValidate
    selectMainField = mainfield{ii};
    for jj = 1:size(subfield,2)
        selectSubField = subfield{ii,jj};
        selectType = typefield{ii,jj};
        selectSize = sizefield{ii,jj};
        
        if ~isempty(selectSubField)
            selectField = objS.(selectMainField).(selectSubField);
            
            % Check the dimension of the selected field.
            objsize = size(selectField);
            if length(objsize) < length(selectSize)
                objsize((length(objsize)+1):length(selectSize)) = 1;
            end
            
            % Check the size of th selected field.
            for kk = 1:length(selectSize)
                if ~ischar(selectSize{kk})
                    valid = valid && (objsize(kk) == selectSize{kk});
                    if ~errorType && ~valid
                        errorType = 1;
                        errorData = ['objS.' selectMainField '.' selectSubField];
                    end
                else
                    if exist(selectSize{kk},'var')
                        valid = valid && (objsize(kk) == eval(selectSize{kk}));
                        if ~errorType && ~valid
                            errorType = 1;
                            errorData = ['objS.' selectMainField '.' selectSubField];
                        end
                    else
                        %assignin('caller',selectSize{kk},objsize(kk));
                        eval([selectSize{kk} '=' num2str(objsize(kk)) ';']);
                    end
                end
            end
            
            % Check the type of the selected field.
            if isa(selectField,'cell')
                if ~isempty(selectField)
                    selectField = selectField{1};
                    valid = valid && (isa(selectField,selectType) || isa(selectField,'sym'));
                    if ~errorType && ~valid
                        errorType = 2;
                        errorData = ['objS.' selectMainField '.' selectSubField];
                    end
                end
            else
                valid = valid && (isa(selectField,selectType) || isa(selectField,'sym'));
                if ~errorType && ~valid
                    errorType = 2;
                    errorData = ['objS.' selectMainField '.' selectSubField];
                end
            end
        end
    end
end

switch errorType
    case 1
        error('spinw:sw_valid:SizeMismatch',['Input argument size mismatch in: ' errorData]);
    case 2
        error('spinw:sw_valid:TypeMismatch',['Input argument type mismatch in: ' errorData]);
end

end