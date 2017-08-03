classdef cif
    % class handling cif data
    %
    % obj = CIF(source)
    %
    % Any field returns a value belonging to that field name. Use
    % cif.('field-name') for field names that are not valid matlab field
    % names. Cif files containing multiple crystal structure, the first one
    % will be retrieved.
    %
    % Input:
    %
    % source        Source of the .cif file, can be a filename, internet
    %               link or a string containing the content of a .cif file.
    %
    % Example:
    %
    % To test loading a .cif file from the internet use:
    %
    % cryst = cif('https://drive.google.com/uc?export=download&id=0BzFs7CQXhehScVFfbkhrZHo1Z1k');
    %
    % See also cif.fieldnames.
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
                    fprintf(['  <a href="matlab:help @cif">cif</a> object with %d fields,\n'...
                        '  stores all data imported from the Crystallographic Information File,\n'...
                        '  source: <a href="matlab:cd %s">%s</a>.\n\n'],numel(obj.cifdat),fileparts(obj.source),obj.source);
                else
                    fprintf(['  <a href="matlab:help @cif">cif</a> object with %d fields,\n'...
                        '  stores all data imported from the Crystallographic Information File,\n'...
                        '  source: %s.\n\n'],numel(obj.cifdat),obj.source);
                end
            else
                fprintf(['  <a href="matlab:help @cif">cif</a> object with %d fields,\n'...
                    '  stores all data imported from the Crystallographic Information File.\n\n'],numel(obj.cifdat));
            end
            
        end
        
    end
    
end
