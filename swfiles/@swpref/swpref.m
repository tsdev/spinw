classdef swpref < dynamicprops
    %SWPREF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Hidden = true, Access=private)
        Name
        Validation
        Value
        Label
    end
    
    properties (Hidden=true, Access = private)
        props = meta.DynamicProperty.empty
    end
    
    methods
        function obj = swpref()
            data = datastruct;
            sPref = obj.get_set_static();
            
            obj.Name = data.Name;
            obj.Validation = data.Validation;
            obj.Value = data.Value;
            obj.Label = data.Label;
            
            if ~isempty(sPref)
                f = fieldnames(sPref);
                for i = 1:length(f)
                    ind = strcmp(f{i},obj.Name);
                    if any(ind)
                        obj.Value{ind} = sPref.(f{i});
                    end
                end
            else
                f = obj.Name;
                for i = 1:length(f)
                    obj.get_set_static(f{i},obj.Value{i});
                    sPref.(f{i}) = obj.Value{i};
                end
            end
            
            for i = 1:length(data.Name)
                obj.props(i) = addprop(obj,obj.Name{i});
                obj.props(i).SetMethod = @(obj, val) set_data(obj,obj.Name{i}, val);
                obj.props(i).GetMethod = @(obj) get_data(obj,obj.Name{i});
                if isfield(sPref,obj.Name{i})
                    obj.(obj.Name{i}) = sPref.(obj.Name{i});
                else
                    obj.(obj.Name{i}) = obj.Value{i};
                end
            end
        end
        
        function out = getpref(obj,name,varargin)
           out = obj.get(name); 
        end
        
        function varargout = get(obj,names)
            if iscell(names)
                j = 1;
                for i = 1:legth(names)
                    if obj.check_names(names{i})
                        varargout{j} = obj.(names{i}); %#ok<AGROW>
                        j = j+1;
                    else
                        error('spref:GetError','There is no field %s in spref',names{i});
                    end
                end
            else
                if obj.check_names(names)
                    varargout{1} = obj.(names);
                else
                    error('spref:GetError','There is no field %s in spref',names);
                end
            end
        end
        
        
        function set(obj,names,values)
            if iscell(names) && iscell(values)
                if length(names) ~= length(values)
                    error('spref:InputError','Names and Values must have the same length')
                end
                for i = 1:length(names)
                    if obj.check_names(names{i})
                        obj.(names{i}) = values{i};
                    else
                        error('spref:SetError','There is no field %s in spref',names{i});
                    end
                end
            else
                if obj.check_names(names)
                    obj.(names) = values;
                else
                    error('spref:SetError','There is no field %s in spref',names);
                end
            end
        end
    end
    
    methods (Hidden=true, Access = private)
        
        function set_data(obj, name, val)
            if ~obj.check_names(name)
                error('spref:GetError','There is no field %s in spref',name);
            end
            
            idx = strcmp(name,obj.Name);
            if ~isempty(obj.Validation{idx})
                checks = obj.Validation{idx};
                for i = 1:length(checks)
                    feval(checks{i},val);
                end
            end
            obj.get_set_static(name, val);
        end
        
        function val = get_data(obj, name)
            if obj.check_names(name)
                val = obj.get_set_static(name);
            else
                error('spref:GetError','There is no field %s in spref',name);
            end
        end
        
        function valid = check_names(obj,name)
            valid = any(strcmp(name,fieldnames(obj)));
        end
    end
    
    methods (Static, Hidden = true)
        function varargout = get_set_static(name,value)
            persistent sPref
            if nargin == 0
                varargout{1} = sPref;
                return
            elseif nargin == 1
                varargout{1} = sPref.(name);
            else
                sPref.(name) = value;
            end
        end
    end
end

