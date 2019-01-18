classdef cif
    % class handling cif data
    %
    % ### Syntax
    %
    % `obj = cif(source)`
    %
    % ### Description
    %
    % `obj = cif(source)` generates a cif class object. The object returns
    % any field value corresponding to an existing field value in the
    % imported .cif file. Use `cif.('field-name')` for field names that are
    % not valid Matlab variable names. If a cif files contains multiple
    % crystal structure, only the first one will be imported.
    %
    % ### Examples
    %
    % To test loading a .cif file from the internet use:
    %
    % cryst = cif('https://goo.gl/Bncwcn');
    %
    % ### Input Arguments
    %
    % `source`
    % : Location of the .cif file. It can be a filename, internet
    %   link or a string containing the content of a .cif file.
    %
    % ### See Also
    %
    % [cif.fieldnames]
    %
    
    properties (Access = private)
        cifdat
        source
        isfile
    end
    
    methods
        function obj = cif(path)
            obj.cifdat = struct('name','','val',[],'type','');
            
            if nargin == 1
                [obj.cifdat, obj.source, obj.isfile] = obj.importcif(path);
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
        
        function fName = fields(obj)
            % returns all the field names of the cif object
            fName = {obj.cifdat(:).name}';
        end
        function disp(obj)
            if ~isempty(obj.source)
                if obj.isfile
                    fprintf(sw_markdown(sprintf(...
                        ['     `cif` object, [cif] class:\n'...
                         '     `Number of fields:` %d\n'...
                         '     `Source:`           [%s]\n'...
                         ],numel(obj.cifdat),obj.source)));
                else
                    fprintf(sw_markdown(sprintf(...
                        ['     `cif` object, [cif] class:\n'...
                         '     `Number of fields:` %d\n'...
                         '     `Source:`           [%s](%s)\n'...
                         ],numel(obj.cifdat),obj.source,obj.source)));
                end
            else
                    fprintf(sw_markdown(sprintf(...
                        ['     `cif` object, [cif] class:\n'...
                         '     `Number of fields:` %d\n'...
                         ],numel(obj.cifdat))));
            end
            
        end
        
    end
    
end
