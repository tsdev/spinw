function obj = initfield(obj)
% SW_INITFIELD(objS) initializes all subfields of the obj structure to
% their initial values.
%

datstruct = datastruct();
mainfield = datstruct.mainfield;
subfield  = datstruct.subfield;
defval    = datstruct.defval;


for ii = 1:length(mainfield)
    for jj = 1:size(subfield,2)
        if ~isempty(subfield{ii,jj})
            if ~isfield(obj,mainfield{ii}) || ~isfield(eval(['obj.' mainfield{ii}]),subfield{ii,jj})
                obj.(mainfield{ii}).(subfield{ii,jj}) = eval(defval{ii,jj});
            end
        end
    end
end

end