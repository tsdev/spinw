function rPref = getpref(prefName, ~)
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
store = 'spinw_global';

% default values
dn = {      'fid'       'pid'             'expert' 'tag'        'nmesh' 'maxmesh' 'npatch' 'fontsize'};
dv = {      1           feature('getpid') 0        'sw_crystal' 2       6         30       12        };

dl = {...
    'file identifier for text output, default value is 1 (Command Window)'...
    'PID value assigned to the running Matlab session, used to reset all pref after restart'...
    'expert mode (1) gives less warnings (not recommended), default value is 0'...
    'defines the tag property of the crystal structure plot figures'...
    'default number of subdivision of the icosahedron for plotting'...
    'maximum number of subdivision of the icosahedron for plotting'...
    'number of edges for patch object'...
    'fontsize for plotting'...
    };

dPref = struct('val',{},'name',{},'label',{});

[dPref(1:numel(dv),1).name] = dn{:};
[dPref(:).label]            = dl{:};
[dPref(:).val]              = dv{:};

% get stored preferences
sPref = getpref(store);

pidNow = feature('getpid');

if ~isempty(sPref) && sPref.pid~=pidNow
    rmpref(store);
    setpref(store,'pid',pidNow);
    sPref = struct('pid',pidNow);
end

prefName = lower(prefName);

if nargin>0
    if strcmp(prefName,'default')
        % return default preference values
        rPref = dPref;
        return
    end
    
    % if a specific value is requested, check if it exist in the default value
    % list
    iPref = find(strcmp(prefName,{dPref(:).name}),1);
    if isempty(iPref)
        error('swpref:getpref:WrongName','The requested SpinW preference does not exists!');
    end
    
    % if a specific value is requested and it exists, return it
    rPref = dPref(iPref);
    
    if isfield(sPref,prefName)
        rPref.val = sPref.(prefName);
    end
    
    if nargin > 1
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

end