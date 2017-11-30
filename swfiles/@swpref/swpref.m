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
        function obj = swpref(opt)
            
            if nargin > 0
                if ischar(opt) && strcmpi(opt,'default')
                    opt = 0;
                else
                    opt = 1;
                end
            else
                opt = 1;
            end
            
            data = datastruct;
            sPref = obj.get_set_static();
            
            obj.Name = data.Name;
            obj.Validation = data.Validation;
            obj.Value = data.Value;
            obj.Label = data.Label;
            
            if opt && ~isempty(sPref)
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
        
        function disp(obj)
            if verLessThan('MATLAB','8.2')
                isTable = false;
            else
                isTable = true;
            end
            
            this_Name  = obj.Name';
            this_Label = obj.Label';
            this_Value = cell(size(this_Label));
            for i = 1:length(this_Name)
               this_Value{i} = obj.(this_Name{i});
            end
                
            varName = {'Name', 'Value', 'Label'};
            var = {this_Name, this_Value, this_Label};
            
            fprintf('Spref object, swpref class:\n')
            if isTable
                disp(table(var{:},'VariableNames',varName))
            else
                for i = 1:length(this_Name)
                    fprintf('%s:\t%s\n',this_Name{i},this_Value{i});
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
        
        runtests(varargin)
        
    end
    
    methods (Static)
        function out = getpref(name,varargin)
            obj = swpref;
            if nargin == 0
                out = obj;
                return
            end
            if ~obj.check_names(name)
                error('spref:GetError','There is no field %s in spref',name);
            end
            if nargin > 1
                out = obj.get(name);
            else
                out.name = name;
                out.label = obj.Label{strcmp(name,obj.Name)};
                out.val = obj.get(name);
            end
        end
        
        function setpref(name,val)
            if nargin < 2
                error('spref:SetInputError','There needs to be a valid name, value pair.');
            end
            obj = swpref;
            if ~obj.check_names(name)
                error('spref:SetError','There is no field %s in spref',name);
            end
            obj.(name) = val;
        end
    end
end

