function RenameFiles(old_string,new_string,filepath,sw)
% 用MATLAB做的一个修改一个文件夹和它的子文件夹名称的程序，共享之，就是把文件名 中的
% 一个字符串用另一个字符串代替，其中的文件夹中各文件遍历方法我觉得是可以借鉴的。
% 
% 是否修改子文件夹中的文件名由第四个开关变量决定，1则修改，0不修改。
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
            RenameFiles(old_string,new_string,filepath,1);  % 调用自己，遍历文件夹
        end;
    end;
end;
%--------------------------------------------------------------------------
end  %End of function RenameFile