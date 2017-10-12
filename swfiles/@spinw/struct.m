function objS = struct(obj)
% converts properties into struct
% 
% ### Syntax
% 
% `objS = struct(obj)`
% 
% ### Description
% 
% `objS = struct(obj)` converts all public properties of `obj` and saves
% them into `objS` struct.
% 
% ### See Also
% 
% [spinw] \| [spinw.copy]
%

objS   = struct;
fNames = fieldnames(obj);
for ii = 1:length(fNames)
    objS.(fNames{ii}) = obj.(fNames{ii});
end

end % struct
