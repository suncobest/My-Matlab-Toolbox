% This script will retify image sequence
if ~exist('drawIM','var'),
    drawIM = input('Display the original and undistorted images or not? ([]=no, other=yes) ','s');
    drawIM = ~isempty(drawIM);
end;

n_frame = 0;
if exist('imgdir','var') && ~isempty(imgdir),
    if exist('imgbase','var') && ~isempty(imgbase) && exist('imgfmt','var') && ~isempty(imgfmt),
        %file name prefix number of original image
        imPreNum = [imgdir '/' imgbase];
        [n_frame, strnum_frame] = check_image_sequence(imPreNum, imgfmt);
        if n_frame ==0,
            fprintf(1,'No images found! Please relocate!\n');
        end;
    end;
end;

while n_frame == 0,
    [imgbase,imgdir]= uigetfile({'*.jpg;*.bmp;*.png;*.gif;*.tif','Image Files (*.jpg,*.bmp,*.png,*.gif,*.tif)';...
        '*.jpg','JPEG Files(*.jpg)'; '*.bmp','BMP Files(*.bmp)'; '*.png','PNG Files(*.png)'; '*.gif','GIF Files(*.gif)';...
        '*.tif','TIFF Files(*.tif)'; '*.*','All Files(*.*)'},'Select the first frame image!');    % Ñ¡ÔñÍ¼Æ¬Â·¾¶
    if imgbase == 0,
        n_frame = 0;
        fprintf(1,'Please relocate images!\n');
    else
        loc_ext = strfind(imgbase,'.');
        imgfmt = imgbase(loc_ext(end)+1:end);   % format of images
        imgbase = input('Basename of image frames (without number nor suffix): ','s');
        if imgdir(end)=='\' || imgdir(end)=='/',
            imgdir(end) = [];
        end;
        imPreNum = [imgdir '/' imgbase];
        [n_frame, strnum_frame] = check_image_sequence(imPreNum, imgfmt);
    end;
end;

% read the 1st image
ima_name = [imPreNum strnum_frame{1} '.' imgfmt];
if imgfmt(1) == 'p',
    if imgfmt(2) == 'p',
        I = double(loadppm(ima_name));
    elseif imgfmt(2) == 'g',
        I = double(loadpgm(ima_name));
    else
        I = double(imread(ima_name));
    end;
else
    if imgfmt(1) == 'r',
        I = readras(ima_name);
    else
        I = double(imread(ima_name));
    end;
end;
[ny,nx,nc] = size(I);

% Pre-compute the necessary indices and blending coefficients to enable quick rectification:
[~,ind_new,ind_1,ind_2,ind_3,ind_4,a1,a2,a3,a4] = rect_index(zeros(ny,nx),eye(3),fc,cc,kc,alpha_c);

for kk = 1:n_frame,
    ima_name = [imPreNum strnum_frame{kk} '.' imgfmt];
    fprintf(1,'Loading image %s...\n',ima_name);
    if imgfmt(1) == 'p',
        if imgfmt(2) == 'p',
            I = double(loadppm(ima_name));
        elseif imgfmt(2) == 'g',
            I = double(loadpgm(ima_name));
        else
            I = double(imread(ima_name));
        end;
    else
        if imgfmt(1) == 'r',
            I = readras(ima_name);
        else
            I = double(imread(ima_name));
        end;
    end;
    Irect = zeros(ny,nx,nc);
    for ii = 1:nc,
        Iii = I(:,:,ii);
        I2ii = Irect(:,:,ii);
        I2ii(ind_new) = a1 .* Iii(ind_1) + a2 .* Iii(ind_2) + a3 .* Iii(ind_3) + a4 .* Iii(ind_4);
        Irect(:,:,ii) = I2ii;
    end;
    if drawIM,
        figure(2);
        image(uint8(I));
        if nc==1,
            colormap(gray(256));
        end;
        drawnow;
        figure(3);
        image(uint8(Irect));
        if nc==1,
            colormap(gray(256));
        end;
        drawnow;
    end;
    
    save_name = [imgdir '/rect_' imgbase strnum_frame{kk} '.' imgfmt];
    fprintf(1,'Saving undistorted image as %s\n', save_name);
    if imgfmt(1) == 'p',
        if imgfmt(2) == 'p',
            saveppm(save_name,uint8(round(Irect)));
        elseif imgfmt(2) == 'g',
            savepgm(save_name,uint8(round(Irect)));
        else
            imwrite(uint8(round(Irect)),save_name,imgfmt);
        end;
    else
        if imgfmt(1) == 'r',
            writeras(save_name,round(Irect),gray(256));
        else
            imwrite(uint8(round(Irect)),save_name,imgfmt);
        end;
    end;
end;
