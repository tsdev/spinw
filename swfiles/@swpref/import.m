function import(obj,location)
% imports swpref object from file
%
% ### Syntax
%
% 
% `obj = import(obj)`
%
% `obj = import(obj,location)`
%
% ### Description
%
% `obj = import(obj)` loads the preferences given in by the file
% `swprefs.json` in the users home folder. It sets the preferences and
% returns a new preference object.
%
% `obj = import(obj,location)` loads the preferences given in by the file 
% specified by `location`, sets the preferences and returns a new
% preference object.
%
% ### See Also
%
% [swpref.export]
%

if nargin == 1
    location = [userpath filesep 'swprefs.json'];
end

if ~exist(location,'file')
   error('spref:ReadError','Can not find the preference file\n%s',location) 
end

f = fopen(location,'r');
c = onCleanup(@() feval(@fclose,f));
[~] = fgetl(f); 

while ~feof(f)
   line = fgetl(f);
   if strcmp(line,'}')
      break 
   end
   name_val = textscan(line,'%s');
   name = name_val{1}{1}(2:end-2);
   value = name_val{1}{2};
   if strcmp(value(end),',')
       value = value(1:end-1);
   end
   if strcmp(value(end),'''')
      % We have a str or fn 
      value = value(2:end-1); % Strip off the '
      if strcmp(value(1),'@')
         value = str2func(value(2:end)); 
      end
   else
      % We have a numeric
      value = str2double(value);
   end
   obj.(name) = value;
end

end