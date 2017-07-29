%%% This script let the user enter the name of the images (base name, numbering scheme,...)
% Modified from data_calib_no_read.m by zpf.

fprintf(1,'\n------------One Camera with Split Views Initialization------------ \n');
n_ima = 0;
while (n_ima==0),
    fprintf(1,'\n');
    calib_name = input('Basename of images (without number nor suffix): ','s');
    format_image = '0';
    while format_image == '0',
        format_image =  input(['Image format: ( ''n''=''png'', ''b''=''bmp'', ''t''=''tif'',' ...
            '''j''=''jpg'',''g''=''jpeg'',\n''p''=''pgm'', ''m''=''ppm'',[]=''r''=''ras'')'],'s');
        if isempty(format_image),
            format_image = 'ras';
        end;
        if lower(format_image(1)) == 'n',
            format_image = 'png';
        elseif lower(format_image(1)) == 'b',
            format_image = 'bmp';
        elseif lower(format_image(1)) == 't',
            format_image = 'tif';
        elseif lower(format_image(1)) == 'j',
            format_image = 'jpg';
        elseif lower(format_image(1)) == 'g',
            format_image = 'jpeg';
        elseif lower(format_image(1)) == 'p',
            format_image = 'pgm';
        elseif lower(format_image(1)) == 'm',
            format_image = 'ppm';
        elseif lower(format_image(1)) == 'r',
            format_image = 'ras';
        else
            disp('Invalid image format');
            format_image = '0'; % Ask for format once again
        end;
    end;
    %  Checks if the images are there in the direcory.
    [n_ima, strnum_cell] = check_image_sequence(calib_name, format_image);
end;

%%% By default, all the images are active for calibration:
active_images = ones(1,n_ima);

% Reading images:
images_read = active_images;
fprintf(1,'\nChecking directory content for the calibration images (no global image loading in memory efficient mode)\n');
fprintf(1,'Found images: ');
nxy_slots = zeros(2,n_ima);      % store nx and ny
for kk = 1:n_ima,
    string_num = strnum_cell{kk};                           % strnum_cellÎª°´Ë³ÐòÅÅÁÐµÄÍ¼Æ¬ÐòºÅ×Ö·ûcell
    ima_name = [calib_name  string_num '.' format_image];
    
    if exist(ima_name,'file')==2,
        fprintf(1,'%d...',kk);
        if mod(kk,20)==0,
            fprintf(1,'\n');
        end;
        if format_image(1) == 'p',
            if format_image(2) == 'p',
                I = double(loadppm(ima_name));
            elseif format_image(2) == 'g',
                I = double(loadpgm(ima_name));
            else
                I = double(imread(ima_name));
            end;
        else
            if format_image(1) == 'r',
                I = readras(ima_name);
            else
                I = double(imread(ima_name));
            end;
        end;
        [ny,nx,~] = size(I);
        nxy_slots(:,kk) = [nx;ny];
    else
        images_read(kk) = 0;
    end;
end;

if ~isequal(nxy_slots, repmat(nxy_slots(:,1),1,n_ima)),
    fprintf(1,'\nERROR: Calibration Images are assumed to have uniform sizes!\n');
    return;
end;

active_images = images_read;
fprintf(1,'\ndone\n');