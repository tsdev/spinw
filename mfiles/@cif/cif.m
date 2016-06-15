classdef cif
    % cif('path')
    %
    % class handling cif data
    % Any field returns a value belonging to that field name. Use
    % cif.('field-name') for field names that are not valid matlab field
    % names. Cif files containing multiple crystal structure, the first one
    % will be retrieved.
    %
    % path      Path to the cif file to be opened.
    %
    % See also cif.fieldnames.
    
    properties (Access = private)
        cifdat
        source
    end
    
    methods
        function obj = cif(path)
            obj.cifdat = struct('name','','val',[],'type','');
            
            if nargin == 1
                [obj.cifdat, obj.source] = obj.importcif(path);
            end
            
        end
        
        function B = subsref(obj, S)
            B = [];
            switch S.type
                case '.'
                    for ii = 1:numel(obj.cifdat)
                        if strcmp(S.subs, obj.cifdat(ii).name)
                            B = obj.cifdat(ii).val;
                            break;
                        end
                    end
                    
                    if strcmp(S.subs, 'cifdat')
                        B = obj.cifdat;
                    end
            end
            
        end
        
        function fName = fieldnames(obj)
            % returns all the field names of the cif object
            fName = {obj.cifdat(:).name}';
        end
        
        function fName = names(obj)
            % returns all the field names of the cif object
            fName = {obj.cifdat(:).name}';
        end
        function disp(obj)
            fprintf(['  <a href="matlab:help @cif">cif</a> object with %d fields,\n'...
                '  stores all data imported from the Crystallographic Information File,\n'...
                '  source file: <a href="matlab:cd %s">%s</a>.\n\n'],numel(obj.cifdat),fileparts(obj.source),obj.source);
        end
        
    end
    
end
