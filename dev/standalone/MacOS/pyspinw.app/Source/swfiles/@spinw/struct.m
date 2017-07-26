function objS = struct(obj)
% extracts all public properties of spinw object into a struct
%
% objS = STRUCT(obj)
%
% See also SPINW, SPINW.COPY.
%

objS   = struct;
fNames = fieldnames(obj);
for ii = 1:length(fNames)
    objS.(fNames{ii}) = obj.(fNames{ii});
end

end % struct
