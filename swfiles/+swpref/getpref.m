function rPref = getpref(prefName, varargin)
% returns SpinW global preferences
% 
% ### Syntax
% 
% `allPref = swpref.getpref`
% 
% `selPref = swpref.getpref(prefName)`
%
% `val = swpref.getpref(prefName,true)`
%
% `rPref = swpref.getpref('default')`
%
% ### Description
% 
% `allPref = swpref.getpref` returns all preference in a struct where each
% field-value pair corresponds to a prefernce name-value pair.
%
% `selPref = swpref.getpref(prefName)` returns a struct that contains the
% value, name and label of the selected preference.
%
% `val = swpref.getpref(prefName,true)` just returns the stored value
% corresponding to `prefName` preference.
%
% `rPref = swpref.getpref('default')` returns the default preference names,
% values and labels of each preferences.
% 
% {{note The preferences are reset after every restart of Matlab, unlike the
% Matlab built-in preferences that are persistent between Matlab sessions.
% If you want certain preferences to keep after closing matlab, define them
% in the `startup.m` file.}}
%
% ### See Also
% 
% [swpref.setpref]
%

if nargin == 0
    % return all current values
    rPref = swpref.pref([]);
    % convert into single struct
    rPref = [{rPref(:).name};{rPref(:).val}];
    rPref = struct(rPref{:});
else
    rPref = swpref.pref(prefName,'get',varargin{:});
end


end