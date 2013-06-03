classdef cif
    % cif('path')
    %
    % class handling cif data
    % Any field returns a value belonging to that field name. Use
    % cif.('field-name') for field names that are not valid matlab field
    % names.
    %
    % path      Path to the cif file to be opened.
    %
    
    properties (Access = private)
        cifdat
    end
    
    methods
        function obj = cif(path)
            obj.cifdat = struct('name','','val',[],'type','');
            
            if nargin == 1
                obj.cifdat = obj.importcif(path);
            end
            
        end
        
        function B = subsref(obj, S)
            B = [];
            switch S.type
                case '.'
                    for ii = 1:numel(obj.cifdat)
                        if strcmp(S.subs, obj.cifdat(ii).name)
                            B = obj.cifdat(ii).val;
                        end
                    end
                    
                    if strcmp(S.subs, 'cifdat')
                        B = obj.cifdat;
                    end
            end
            
        end
    end
    
end
