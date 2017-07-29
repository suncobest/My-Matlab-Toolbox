% read txt

file=fopen('Grid_Data.txt','r');
data=[];
while 1
    line=fgetl(file);
    if ~ischar(line)
        break;
    end
    line=str2num(line);
    data=[data;line];
end
fclose(file);