function writex3dfile(folder, filename, data)
% Convert the XML data to a string and write it to a X3D file
[~,~,str]=XMLtostring(data);
fid = fopen([folder filename '.x3d'],'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<!DOCTYPE X3D PUBLIC "ISO//Web3D//DTD X3D 3.1//EN" "http://www.web3d.org/specifications/x3d-3.1.dtd">\n');
fprintf(fid,str);
fclose(fid);

