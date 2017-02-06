function rPref = pref(prefName, mode, value)
% returns SpinW global preferences
%
% rPref = swpref.getpref
%
% The preferences are reset after every restart of Matlab, unlike the
% Matlab built-in preferences that are persistent between Matlab sessions.
% If you want certain preferences to keep after closing matlab, define them
% in the <a href="matlab:edit('startup.m')">startup.m</a> file.
%
% swpref.getpref() returns the names, values and labels of each
% preferences. Default values are returned, if no values are saved. rPref
% is a struct with field names 'name', 'label' and 'val'. Each field is a
% cell.
%
% rPref = swpref.getpref(pName, {simple})
%
% Returns only the requested SpinW preference name, value and label in a
% struct. Each field contains the requested value. If a second argument is
% given (simple) with any value, only the value of the preference is
% returned.
%
% rPref = swpref.getpref('default')
%
% Returns the default names, values and labels of each preferences.
%
% See also GETPREF, SETPREF, SWPREF.SETPREF.

% the storage name within built-in getpref/setpref
%store = 'spinw_global';

prefName = lower(prefName);

% default values
dn = {'fid' 'expert' 'tag'    'nmesh' 'maxmesh' 'npatch' 'fontsize' 'tid' 'colormap' };
dv = {1     0        'swplot' 1       6         20       12         1     @cm_inferno};

dl = {...
    'file identifier for text output, default value is 1 (Command Window)'...
    'expert mode (1) gives less warnings (not recommended), default value is 0'...
    'defines the tag property of the crystal structure plot figures'...
    'default number of subdivision of the icosahedron for plotting'...
    'maximum number of subdivision of the icosahedron for plotting'...
    'number of edges for patch object'...
    'fontsize for plotting'...
    'identifier how the timer is shown, default value is 1 (Command Window), value 2 gives graphical output'...
    'default colormap'...
    };

dPref = struct('val',{},'name',{},'label',{});

[dPref(1:numel(dv),1).name] = dn{:};
[dPref(:).label]            = dl{:};
[dPref(:).val]              = dv{:};

% store the preferences in a persistent variable
persistent sPref
% lock the file in memory
mlock

if isempty(sPref)
    % start with the default values
    sPref = dPref;
end

if isempty(prefName)
    % just return all current values
    rPref = sPref;
    return
end

switch mode
    case 'get'
                
        if ~isempty(prefName)
            if strcmp(prefName,'default')
                % return default preference values
                rPref = dPref;
                return
            end
            
            % if a specific value is requested, check if it exist in the default value
            % list
            iPref = find(strcmp(prefName,{dPref(:).name}),1);
            if isempty(iPref)
                error('pref:WrongName','The requested SpinW preference does not exists!');
            end
            
            % if a specific value is requested and it exists, return it
            rPref = sPref(iPref);
            
            if nargin>2
                rPref = rPref.val;
            end
            
            return
        else
            % return all stored values
            rPref = dPref;
            % overwrite default values for existing preferences
            if ~isempty(sPref)
                fPref = fieldnames(sPref);
                for ii = 1:numel(fPref)
                    rPref(strcmp(fPref{ii},{dPref(:).name})).val = sPref.(fPref{ii});
                end
            end
            
        end
    case 'set'
        % set preferences
        if strcmp(prefName,'default')
            sPref = dPref;
            return
        end
                
        % check if the preference name exists
        iPref = find(strcmp(prefName,{dPref(:).name}),1,'first');
        if ~isempty(iPref)
            % check if the preferences label contains a choice string
            %str0 = strsplit(dPref(iPref).label,' ');
            str0 = regexp(dPref(iPref).label,' ','split');
            %opt0 = strsplit(str0{end},'/');
            opt0 = regexp(str0{end},'/','split');
            
            if numel(opt0) > 1
                % there is a choice of different string options
                if ~ischar(value) || ~any(strcmp(value,opt0))
                    error('pref:WrongInput',['The selected preference has a restricted choice: ' str0{end} '!'])
                end
                sPref(iPref).val = value;
            else
                % the value has to be a scalar
                % TODO check for other type of values
                sPref(iPref).val = value;
            end
            
        else
            error('pref:WrongName','The given name is not a valid SpinW global preferences!');
        end
end

end