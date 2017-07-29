%%% This script let the user enter the name of the images (base name, numbering scheme,...)
% Modified from data_calib_no_read.m by zpf.

fprintf(1,'\n----------------Multicams  Initialization---------------- \n');
if ~exist('n_cam','var'),
    n_cam = input('Number of cameras = ([] = 1)');
    if isempty(n_cam),
        n_cam = 1;
    end;
end;
fprintf(1,'\nThere should be identical number of images in every camera view!\n');

% hand_list(i) means axis handness of the i-th camera view
if ~exist('hand_list','var'),
    flag = 1;
    fprintf(1,'\nHandedness vector of cameras do not exist!\n');
else
    hand_list =  hand_list(:)';
    flag = ~isequal(abs(hand_list),ones(1,n_cam));
    if flag,
        fprintf(1,'\nUnexpected value for handedness vector! Please input again!\n');
    end;
end;
while flag,
    fprintf(1,['\nPlease input the axis handedness vector of all %d camera views!\n' ...
        'Note: 1 stands for right-handed, -1 stands for left-handed.\n'],n_cam);
    hand_list = input(['handedness vector = ([] = [' num2str(ones(1,n_cam)) '])']);
    if isempty(hand_list),
        hand_list = ones(1,n_cam);
        flag = 0;
    else
        hand_list =  hand_list(:)';
        flag = ~isequal(abs(hand_list),ones(1,n_cam));
        if flag,
            fprintf(1,'\nUnexpected value for handedness vector! Please input again!\n');
        end;
    end;
end;

nima_slot = zeros(1,n_cam);
imbase = cell(1,n_cam);
imformat = imbase;
imstrnum = imbase;
diff_imnum = imbase;      % difference of image numbers
for pp = 1:n_cam,
    n_ima = 0;
    while n_ima==0,
        fprintf(1,'\n');
        calib_name = input(['Image basename of camera ' num2str(pp) ' (without number nor suffix): '],'s');
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
        [n_ima, strnum_cell, image_num] = check_image_sequence(calib_name, format_image);
    end;
    nima_slot(pp) = n_ima;
    diff_imnum{pp} = diff(image_num);
    imbase{pp} = calib_name;
    imformat{pp} = format_image;
    imstrnum{pp} = strnum_cell;
end;

if ~all(nima_slot==n_ima),
    fprintf(1,'\nERROR: Image numbers of every camera view are not consistent with each other!\n');
    return;
end;
if n_ima>1,
    for pp = 2:n_cam,
        if ~isequal(diff_imnum{pp}, diff_imnum{1}),
            fprintf(1,'\nERROR: Image number differences of camera %d are not consistent with camera 1!\n',pp);
            return;
        end;
    end;
end;

%%% By default, all the images are active for calibration:
active_images = ones(1,n_ima);
% active_imgviews(i,j) denote the availability of the j-th view of i-th image:1 for yes£¬0 for no
active_imgviews = ones(n_cam,1) * active_images;

% Reading images:
imsize = zeros(2,n_cam);
images_read = active_images;
fprintf(1,'\nChecking directory content for the calibration images (memory efficient mode)\n');
for pp = 1:n_cam,
    calib_name = imbase{pp};
    format_image = imformat{pp};
    strnum_cell = imstrnum{pp};
    nxy_slots = zeros(2,n_ima);      % store nx and ny
    fprintf(1,'Found images of Camera %d: ',pp);
    for kk = 1:n_ima,
        string_num = strnum_cell{kk};    % strnum_cellÎª°´Ë³ÐòÅÅÁÐµÄÍ¼Æ¬ÐòºÅ×Ö·ûcell
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
        fprintf(1,'\nERROR: Images of camera %d have different sizes!\n', pp);
        return;
    end;
    imsize(:,pp) = [nx;ny];
    fprintf(1,['\nDone with Cam' num2str(pp) '.\n']);
    active_imgviews(pp,:) = images_read;
end;

active_images = any(active_imgviews,1);
ind_active = find(active_images);
fprintf(1,'\nAll done!\n');
