function str=base64_encode(filename)
% str = base64_encode(filename)
%

% Base64 Characters
dstr='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';

% Open File
fid = fopen(filename,'r', 'ieee-be');

% Determine length of file
fseek(fid,0,'eof');
fsize=ftell(fid);
fseek(fid,0,'bof');

% Get all chunks of 64bit (4*6bit = 3 bytes = 3*8 bit)
fr=4*floor(fsize/3);
c = fread(fid, fr, 'ubit6=>uint8');
% Get the remaining 1 or 2 bytes (8bit or 16bit)
d = fread(fid, inf, 'ubit1=>uint8');
fclose(fid);

% Determine Preamble, and convert also the last bits to 6bits by adding
% a few extra zero bits
if(~isempty(d))
    switch length(d)
        case 8
            d(9:12)=0;
            d=reshape(d,[6 2]);
            pr='==';
        case 16
            d(17:18)=0;
            d=reshape(d,[6 3]);
            pr='=';
    end
    d(1,:)=d(1,:)*32; d(2,:)=d(2,:)*16; d(3,:)=d(3,:)*8;
    d(4,:)=d(4,:)*4;  d(5,:)=d(5,:)*2; d(6,:)=d(6,:)*1;
    d=sum(d,1)';
    c=[c;d];
else
    pr='';
end

% Convert the 64bit values to Base64
str=[dstr(c+1) pr];

