function setpref(prefName, varargin)
% sets SpinW global preferences
% 
% ### Syntax
% 
% `swpref.setpref(prefName, value)`
% 
% `swpref.setpref('default')`
%
% ### Description
% 
% `swpref.setpref(prefName, value)` sets the value of `prefName`
% preferences.
%
% `swpref.setpref('default')` resets all preference values to the default one.
%
% {{note The preferences are reset after every restart of Matlab, unlike the
% Matlab built-in preferences that are persistent between Matlab sessions.
% If you want certain preferences to keep after closing matlab, define them
% in the `startup.m` file.}}
% 
% ### See Also
% 
% [swpref.getpref]
%

if nargin>0 && mod(nargin,2)==0
    % check for usemex option
    mexopt = find(strcmpi('usemex',varargin));
    
    if ~isempty(mexopt) && varargin{mexopt+1}
        % check for the existence of the necessary mex files
        if exist('chol_omp','file')==3 && exist('eig_omp','file')==3
        else
            warning('setpref:MissingMex','Necessary mex files are missing, compile them!')
            varargin{mexopt+1} = false;
        end
    end
end

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