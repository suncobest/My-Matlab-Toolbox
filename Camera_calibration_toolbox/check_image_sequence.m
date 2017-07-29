function [n_ima, strnum_cell, image_num, ind_valid] = check_image_sequence(base_name, format_image)

% Check image sequence in the current direcory and.
% Sort image numbers and store them in image_num and strnum_cell (string format)
% Important value: n_ima,strnum_cell,image_num,ind_valid

imagedir = dir([base_name '*.' format_image]);
Nl = size(imagedir,1);
n_ima = 0;
ind_valid = [];
if Nl == 0,
    fprintf(1,'No image found. Basename may be wrong.\n');
    return;
end;

loc_base = find(base_name=='\' | base_name=='/',1,'last');
if ~isempty(loc_base),
    base_name = base_name(loc_base(end)+1 : end);
end;

loc_extension = [];
length_name = size(base_name,2);
for kk = 1:Nl,
    filenamepp =  imagedir(kk).name;
    if isempty(base_name),
        loc_base = 1;
    else
         % strfind(TEXT,PATTERN)返回字符串TEXT中PATTERN出现的第一个指标，出现几次返回几个
        loc_base = strfind(filenamepp,base_name);
    end;
    % loc_base取第1个字符，loc_ext是文件扩展名的首字符位置，取最后1个
    loc_ext = strfind(filenamepp,format_image);
    % string_num为通用名calib_name和扩展名format_image之间的数字（去除'.'）
    string_num = filenamepp(length_name+1:loc_ext(end) - 2);
    
    if isempty(str2double(string_num)),
        loc_base = [];
    end;
    if ~isempty(loc_base),
        if (loc_base(1) ~= 1),
            loc_base = [];
        end;
    end;
    if ~isempty(loc_base) && ~isempty(loc_ext),
        n_ima = n_ima + 1;
        ind_valid = [ind_valid kk];                             % 合乎要求的图片在imagedir中的序号
        loc_extension = [loc_extension loc_ext(end)];           % 合乎要求的图片扩展名首字符的位置列表
    end;
end;

if n_ima==0,
    % try the upper case format:
    format_image = upper(format_image);
    imagedir = dir([base_name '*.' format_image]);
    for kk = 1:Nl,
        filenamepp =  imagedir(kk).name;
        if isempty(base_name),
            loc_base = 1;
        else
            loc_base = strfind(filenamepp,base_name);
        end;
        loc_ext = strfind(filenamepp,format_image);
        string_num = filenamepp(length_name+1:loc_ext(end) - 2);
        
        if isempty(str2double(string_num)),
            loc_base = [];
        end;
        if ~isempty(loc_base),
            if (loc_base(1) ~= 1),
                loc_base = [];
            end;
        end;
        if ~isempty(loc_base) && ~isempty(loc_ext),
            n_ima = n_ima + 1;
            ind_valid = [ind_valid kk];
            loc_extension = [loc_extension loc_ext(end)];
        end;
    end;
    
    if (n_ima==0),
        fprintf(1,'No image found. File name may be wrong.\n');
        return;
    end;
end;

% Get all the string numbers:
strnum_cell = cell(1,n_ima);                         % strnum_cell为cell型变量，存储图片的序号字符
image_num =  zeros(1,n_ima);

for kk = 1:n_ima,
    name = imagedir(ind_valid(kk)).name;
    string_num = name(length_name+1:loc_extension(kk) - 2);
    strnum_cell{kk} = string_num;
    image_num(kk) = str2double(string_num);     % image_num为图片的序数列表
end;

% sort image number
[image_num,ind] = sort(image_num);                               % 将图片序号由小到大排列
strnum_cell = strnum_cell(ind);                          % 得到按小到大排列的序号字符列表cell
ind_valid = ind_valid(ind);                              % 得到按小到大序号排列的有效图片在imagedir中的指标

return;
% 得到按小到大序号排列的有效图片的扩展名首字符位置列表
% loc_extension = loc_extension(ind);           
