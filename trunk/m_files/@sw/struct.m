function objS = struct(obj)
% extracts all public properties of sw object into a struct
%

objS   = struct;
fNames = fieldnames(obj);
for ii = 1:length(fNames)
    objS.(fNames{ii}) = obj.(fNames{ii});
end

end % struct
