function setpref(prefName, varargin)
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

if nargin<=2
    swpref.pref(prefName,'set',varargin{:});
elseif mod(nargin,2)==0
    arg = [{prefName} varargin];
    % multiple variable
    for ii = 1:2:numel(arg)
        swpref.pref(arg{ii},'set',arg{ii+1});
    end
else
    error('setpref:WrongNumberOfArgument','Wrong number of input arguments!')
end

end