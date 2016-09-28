function setpref(prefName, value)
% sets SpinW global preferences
%
% swpref.setpref(prefName, value)
%
% The preferences are reset after every restart of Matlab, unlike the
% Matlab built-in preferences that are persistent between Matlab sessions.
% If you want certain preferences to keep after closing matlab, define them
% in the <a href="matlab:edit('startup.m')">startup.m</a> file.
%
% swpref.setpref() sets the value of the prefName in the SpinW global
% preferences.
%
% swpref.setpref('default')
%
% Resets all preference values to the default one.
%
% See also SWPREF.SETPREF.
%

% the storage name within built-in getpref/setpref
store = 'spinw_global';

pidNow = feature('getpid');

if ispref(store) && getpref(store,'pid')~=pidNow
    rmpref(store);
end

setpref(store,'pid',pidNow);

if strcmp(prefName,'default')
    if ispref(store)
        rmpref(store);
    end
    setpref(store,'pid',pidNow);
    return
end

if strcmp(prefName,'pid')
    warning('swpref:setpref:Locked','pid value can not be changed!')
    return
end

% check if the preference name exists
dPref = swpref.getpref('default');

iPref = find(strcmp(prefName,{dPref(:).name}),1,'first');
if ~isempty(iPref)
    % check if the preferences label contains a choice string
    str0 = strsplit(dPref(iPref).label,' ');
    opt0 = strsplit(str0{end},'/');
    
    if numel(opt0) > 1
        % there is a choice of different string options
        if ~ischar(value) || ~any(strcmp(value,opt0))
            error('swpref:setpref:WrongInput',['The selected preference has a restricted choice: ' str0{end} '!'])
        end
        setpref(store,prefName,value);
    else
        % the value has to be a scalar
        % TODO check for other type of values
        setpref(store,prefName,value);
    end
    
else
    error('swpref:setpref:WrongName','The given name is not a valid SpinW global preferences!');
end

end