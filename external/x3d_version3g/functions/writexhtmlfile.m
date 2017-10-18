function writexhtmlfile(folder, filename, data)
% Convert the XML data to a string and write it to a XHTML file
[~,~,str]=XMLtostring(data);
fid = fopen([folder filename '.xhtml'],'w');
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\n');
fprintf(fid,str);
fclose(fid);

% Create the supporting X3DOM  files
%load('x3dom.mat')
%folder=[folder 'x3dom'];
%if(~isdir(folder)), mkdir(folder); end
%for i=1:length(data)
%    fid = fopen([folder '/' data(i).filename], 'w','ieee-le');
%    fwrite(fid,data(i).filedata,'uint8');
%    fclose(fid);
%end