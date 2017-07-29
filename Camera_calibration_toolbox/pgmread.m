function img = pgmread(filename)
% function img = pgmread(filename)
%   this is my version of pgmread for the pgm file created by XV.
%  this program also corrects for the shifts in the image from pm file.

fid = fopen(filename,'r');
fscanf(fid, 'P5\n');
cmt = '#';
while strfind(cmt, '#'),
  cmt = fgets(fid);
   if length(strfind(cmt, '#')) ~= 1,
      YX = sscanf(cmt, '%d %d');
      y = YX(1); x = YX(2);
   end;
end;

img = fread(fid,[y,x],'uint8');
img = img';
fclose(fid);

