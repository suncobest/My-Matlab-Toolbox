function RenameFiles(old_string,new_string,filepath,sw)
% ��MATLAB����һ���޸�һ���ļ��к��������ļ������Ƶĳ��򣬹���֮�����ǰ��ļ��� �е�
% һ���ַ�������һ���ַ������棬���е��ļ����и��ļ����������Ҿ����ǿ��Խ���ġ�
% 
% �Ƿ��޸����ļ����е��ļ����ɵ��ĸ����ر���������1���޸ģ�0���޸ġ�
%Rename files

%--------------------------------------------------------------------------
%Rename subfolder files or not
if nargin<4 || isempty(sw),
    sw=0;
else
    sw=~~sw;
end;

if nargin<3 || isempty(filepath),
    filepath = pwd;
end;
%filepath initialization   
if filepath(end)=='\' || filepath(end)=='/',
    filepath=filepath(1:end-1);
end;
%--------------------------------------------------------------------------
%Rename current folder
files=dir(filepath);
if ispc,
    for i=1:length(files),
        filename=files(i).name;
        if ~isempty(strfind(filename,old_string)),
            newfilename=strrep(filename,old_string,new_string);
            dos(['rename "' filepath, '\', filename, '" "', newfilename, '"']);
        end;
    end;
else
    for i=1:length(files),
        filename=files(i).name;
        if ~isempty(strfind(filename,old_string)),
            newfilename=strrep(filename,old_string,new_string);
            movefile([filepath,'/',filename],[filepath,'/',newfilename]);
        end;
    end;
end;
%--------------------------------------------------------------------------
%Rename subfolder
if sw,
    for i=1:length(files),
        filename=files(i).name;
        if isdir(filename)==1 && filename~='.' && ~strcmp(filename,'..'),
            filepath=cat(2,filepath,'/',filename);
            RenameFiles(old_string,new_string,filepath,1);  % �����Լ��������ļ���
        end;
    end;
end;
%--------------------------------------------------------------------------
end  %End of function RenameFile