classdef swpref < dynamicprops
    % class to store and retrieve persistent settings
    %
    % ### Syntax
    %
    % `pref = swpref`
    %
    % `pref = swpref('default')`
    %
    % ### Description
    %
    % `pref = swpref` retrieves and creates a preference object.
    %
    % `pref = swpref('default')` resets all preferences to default values.
    %
    % The settings sotred in the `swpref` class for spinw objects will
    % persist during a single Matlab session. It is different from the
    % Matlab built-in preferences, as swpref resets all settings to factory
    % default after every restart of Matlab.
    %
    % ### Examples
    %
    % We change the fontsize value and show that it is retained even when a
    % new instance of the object is created:
    %
    % ```
    % >>pref = swpref
    % >>pref.fontsize>>
    % >>pref.fontsize = 18
    % >>pref2 = swpref
    % >>pref.fontsize>>
    % >>pref2.fontsize>>
    % ```
    %
    % ### Properties
    %
    % Properties can be changed by directly assigning a new value to them.
    % Once a new value to a given property is assigned, it will be retained
    % until the end of a MATLAB session, even if a new class instance is
    % created.
    %
    % ### Methods
    %
    % Methods are the different commands that require an `swpref` object as
    % a first input, thus they can be called as `method1(obj,...)`,
    % alternatively the equivalent command is `obj.method1(...)`.
    %
    % swpref.get
    % swpref.set
    % swpref.export
    % swpref.import
    %
    % ### Commands
    %
    % Commands are methods which can be called without first creating a
    % preference object 'swpref.command(....)'.
    %
    % swpref.getpref
    % swpref.setpref
    %
    
    properties(Hidden = true, Access=private)
        % stores the details to create and check dynamic properties.
        %
        % `Name` a cell array giving the name of each dynamic property.
        %
        % `Validation` a cell array with functions to evaluate when a
        % property is set.
        %
        % `Value` a cell array giving the default value of a dynamic property
        %
        % `Label` a cell array giving a short description on the dynamic
        % property.
        %
        % These details are retrieved from the private file `datastruct.m`.
        %
        Name
        Validation
        Value
        Label
    end
    
    properties (Hidden=true, Access = private)
        % stores the details of the dynamic properties.
        props = meta.DynamicProperty.empty
    end
    
    methods
        function obj = swpref(opt)
            % Spin preference constructor.
            %
            % ### Syntax
            %
            % `pref = swpref`
            %
            % `pref = swpref('default')`
            %
            %
            % ### Description
            %
            % `pref = swpref` retrieves and creates a preference object.
            %
            % `pref = swpref('default')` resets all preferences to default values.
            %
            %
            % {{note The preferences are reset after every restart of Matlab, unlike the
            % Matlab built-in preferences that are persistent between Matlab sessions.
            % If you want certain preferences to keep after closing matlab, use the
            % 'pref.export(fileLocation)' and 'pref.import(fileLocation)' functions.
            % These can be added to your startup file.}}
            %
            % ### See Also
            %
            % [swpref.get], [swpref.set]
            %
            
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
            % retrieves a preference value
            %
            % ### Syntax
            %
            % `value = get(obj, name)`
            %
            % `value = obj.get(name)`
            %
            % ### Description
            %
            % `value = get(obj, name)` gets the preference `name`.
            %
            % `value = obj.get(name)` gets the preference `name`.
            %
            % ### See Also
            %
            % [swpref.set]
            %
            
            if nargin == 1
                error('swpref:GetError','You need to supply a parameter to get!');
            end
            
            if iscell(names)
                j = 1;
                for i = 1:legth(names)
                    if obj.check_names(names{i})
                        varargout{j} = obj.(names{i}); %#ok<AGROW>
                        j = j+1;
                    else
                        error('swpref:GetError','There is no field %s in swpref',names{i});
                    end
                end
            else
                if obj.check_names(names)
                    varargout{1} = obj.(names);
                else
                    error('swpref:GetError','There is no field %s in swpref',names);
                end
            end
        end
        
        
        function set(obj,names,values)
            % sets a preference value
            %
            % ### Syntax
            %
            % `set(obj, name, value)`
            %
            % `obj.set(name, value)`
            %
            % ### Description
            %
            % `set(obj, name, value)` sets the preference `name` to the
            % value given by `value`
            %
            % `obj.set(name, value)` sets the preference `name` to the
            % value given by `value`
            %
            % ### See Also
            %
            % [swpref.get]
            %
            
            if nargin < 2
                error('swpref:SetError','You need to supply a parameter,value pair to set!');
            end
            
            if iscell(names) && iscell(values)
                if length(names) ~= length(values)
                    error('swpref:SetInputError','Names and Values must have the same length')
                end
                for i = 1:length(names)
                    if obj.check_names(names{i})
                        obj.(names{i}) = values{i};
                    else
                        error('swpref:SetError','There is no field %s in swpref',names{i});
                    end
                end
            else
                if obj.check_names(names)
                    obj.(names) = values;
                else
                    error('swpref:SetError','There is no field %s in swpref',names);
                end
            end
        end
        
        function disp(obj)
            % function called when a preference is displayed.
            %
            % ### Syntax
            %
            % `disp(obj)`
            %
            % `obj.disp`
            %
            % ### Description
            %
            % `disp(obj)` shows a table of all dynamic properties including
            % name, value and description. If a table is not supported only
            % `name: value` is displayed for all properties.
            %
            
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
            
            fprintf(sw_markdown(sprintf(...
                ['     `Swpref` object, [swpref] class:\n' ...
                '     `Stored preferences:`\n'])))
            if isTable
                disp(table(var{:},'VariableNames',varName))
            else
                disp(swpref.getpref);
            end
        end
    end
    
    methods (Hidden=true, Access = private)
        
        function set_data(obj, name, val)
            % Function called when a vairable is set.
            %
            %  {{warning Internal function for the Spin preferences.}}
            %
            % ### Syntax
            %
            % 'set_data(obj, name, value)'
            %
            % ### Description
            %
            % 'set_data(obj, name, value)' sets the 'value' of parameter
            % 'name' which is stored in persistent storage.
            %
            % ### See Also
            %
            % [swpref.setpref], [swpref.get_data]
            %
            
            if ~obj.check_names(name)
                error('swpref:SetError','There is no field %s in swpref',name);
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
            % Function called when a vairable is retrieved.
            %
            %  {{warning Internal function for the Spin preferences.}}
            %
            % ### Syntax
            %
            % 'value = get_data(obj, name)'
            %
            % ### Description
            %
            % 'value = get_data(obj, name)' returns the value of parameter
            % 'name' from persistent storage.
            %
            % ### See Also
            %
            % [swpref.setpref], [swpref.set_data]
            %
            
            if obj.check_names(name)
                val = obj.get_set_static(name);
            else
                error('swpref:GetError','There is no field %s in swpref',name);
            end
        end
        
        function valid = check_names(obj,name)
            % Checking to see if a get/set name is valid.
            %
            %  {{warning Internal function for the Spin preferences.}}
            %
            % ### Syntax
            %
            % 'logical = obj.check_names(name)'
            %
            % ### Description
            %
            % 'logical = obj.check_names(name)' returns true if 'name' is a
            % valid field of 'obj' and false otherwise.
            %
            
            valid = any(strcmp(name,fieldnames(obj)));
        end
    end
    
    methods (Static, Hidden = true)
        function varargout = get_set_static(name,value)
            % The internal session persistent storage of vairables.
            %
            %  {{warning Internal function for the Spin preferences.}}
            %
            % ### Syntax
            %
            % 'values = get_set_static()'
            %
            % 'value' = get_set_static(prefName)'
            %
            % 'get_set_static(prefName, value)'
            %
            % ### Description
            %
            % 'values = get_set_static()' retrieves all the preferences in
            % the storage.
            %
            % 'value' = get_set_static(prefName)' returns the value of
            % preference given by 'prefName'.
            %
            % 'get_set_static(prefName, value)' sets the preference given
            % by 'prefName' to value 'value'
            %
            % ### See Also
            %
            % [swpref.check_names], [swpref.get_data], [swpref.set_data]
            %
            
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
    
    methods (Static)
        function out = getpref(name,varargin)
            % returns SpinW global preferences
            %
            % ### Syntax
            %
            % `allPref = swpref.getpref`
            %
            % `selPref = swpref.getpref(prefName)`
            %
            % `val = swpref.getpref(prefName,[])`
            %
            % ### Description
            %
            % `allPref = swpref.getpref` returns all preference in a struct where each
            % field-value pair corresponds to a prefernce name-value pair.
            %
            % `selPref = swpref.getpref(prefName)` returns a struct that contains the
            % value, name and label of the selected preference.
            %
            % `val = swpref.getpref(prefName,[])` just returns the stored value
            % corresponding to `prefName` preference.
            %
            % {{note The preferences are reset after every restart of Matlab, unlike the
            % Matlab built-in preferences that are persistent between Matlab sessions.
            % If you want certain preferences to keep after closing matlab, define them
            % in the `startup.m` file.}}
            %
            % {{warning This is a legacy function. It is better to use `pref = swpref` and then
            % use regular `pref.prefName = value` or `set(pref, prefName, value)` syntax.}}
            %
            % ### See Also
            %
            % [swpref.setpref]
            %
            
            obj = swpref;
            if nargin == 0
                val = cellfun(@(C)obj.(C),obj.Name,'UniformOutput',false);
                out = [obj.Name; val];
                out = struct(out{:});
                return
            end
            if ~obj.check_names(name)
                error('swpref:GetError','There is no field %s in swpref',name);
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
            % sets SpinW global preferences
            %
            % ### Syntax
            %
            % `swpref.setpref(prefName, value)`
            %
            % ### Description
            %
            % `swpref.setpref(prefName, value)` sets the value of `prefName`
            % preferences.
            %
            % {{note The preferences are reset after every restart of Matlab, unlike the
            % Matlab built-in preferences that are persistent between Matlab sessions.
            % If you want certain preferences to keep after closing matlab, define them
            % in the `startup.m` file.}}
            %
            % {{warning This is a legacy function. It is better to use 'pref = swpref' and then
            % use regular 'value = pref.prefName' or 'get(pref, prefName)' syntax.}}
            %
            % ### See Also
            %
            % [swpref.getpref]
            %
            
            if nargin < 2
                error('swpref:SetInputError','There needs to be a valid name, value pair.');
            end
            obj = swpref;
            if ~obj.check_names(name)
                error('swpref:SetError','There is no field %s in swpref',name);
            end
            obj.(name) = val;
        end
        
        function export(varargin)
            % saves swpref object into a file
            %
            % ### Syntax
            %
            % `success = swpref.export(location)`
            %
            % `success = swpref.export`
            %
            % ### Description
            %
            % `success = swpref.export(location)` writes the preferences given by swpref to
            % a file location given by `location`. The file is in a basic `.json`
            % format.
            %
            % `success = swpref.export` writes the preferences given in `obj` to
            % the users home folder as `swprefs.json`. The file is in a basic `.json`
            % format.
            %
            % ### See Also
            %
            % [swpref.import]
            %
            obj = swpref;
            obj.export(varargin{:})
        end
        
        
        function obj = import(varargin)
            % imports swpref object from file
            %
            % ### Syntax
            %
            %
            % `obj = swpref.import`
            %
            % `obj = swpref.import(location)`
            %
            % ### Description
            %
            % `obj = swpref.import` loads the preferences given in by the file
            % `swprefs.json` in the users home folder. It sets the preferences and
            % returns a new preference object.
            %
            % `obj = swpref.import(location)` loads the preferences given in by the file
            % specified by `location`, sets the preferences and returns a new
            % preference object.
            %
            % ### See Also
            %
            % [swpref.export]
            %
            
            obj = swpref;
            obj = obj.import(varargin{:});
        end
        
    end
end
