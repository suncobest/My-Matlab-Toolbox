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
         % strfind(TEXT,PATTERN)�����ַ���TEXT��PATTERN���ֵĵ�һ��ָ�꣬���ּ��η��ؼ���
        loc_base = strfind(filenamepp,base_name);
    end;
    % loc_baseȡ��1���ַ���loc_ext���ļ���չ�������ַ�λ�ã�ȡ���1��
    loc_ext = strfind(filenamepp,format_image);
    % string_numΪͨ����calib_name����չ��format_image֮������֣�ȥ��'.'��
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
        ind_valid = [ind_valid kk];                             % �Ϻ�Ҫ���ͼƬ��imagedir�е����
        loc_extension = [loc_extension loc_ext(end)];           % �Ϻ�Ҫ���ͼƬ��չ�����ַ���λ���б�
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
strnum_cell = cell(1,n_ima);                         % strnum_cellΪcell�ͱ������洢ͼƬ������ַ�
image_num =  zeros(1,n_ima);

for kk = 1:n_ima,
    name = imagedir(ind_valid(kk)).name;
    string_num = name(length_name+1:loc_extension(kk) - 2);
    strnum_cell{kk} = string_num;
    image_num(kk) = str2double(string_num);     % image_numΪͼƬ�������б�
end;

% sort image number
[image_num,ind] = sort(image_num);                               % ��ͼƬ�����С��������
strnum_cell = strnum_cell(ind);                          % �õ���С�������е�����ַ��б�cell
ind_valid = ind_valid(ind);                              % �õ���С����������е���ЧͼƬ��imagedir�е�ָ��

return;
% �õ���С����������е���ЧͼƬ����չ�����ַ�λ���б�
% loc_extension = loc_extension(ind);           
