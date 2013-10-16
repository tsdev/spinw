function varargout = sw_version()
% SW_VERSION gives the current version of the SpinW code
%

% read file header
testStr = '% $Author: Sandor Toth$ ($Contact: sandor.toth@psi.ch$)';

while numel(testStr) > 0
    [~, testStr] = strtok(testStr,'$');
    [partStr, testStr] = strtok(testStr,'$');
    partStr
end



build   = 'SpinW release 2.0beta';
name    = 'SpinW';
version = 2.0;

build_data = '15.10.2013';
contact    = 'sandor.toth@psi.ch';
licence    = 'GNU GENERAL PUBLIC LICENSE';

if nargout ~= 1
    disp(build)
else
    varargout{1}.name = name;
    varargout{1}.version = version;
    varargout{1}.build_data = build_data;
    varargout{1}.contact = contact;
    varargout{1}.licence = licence;
end

end