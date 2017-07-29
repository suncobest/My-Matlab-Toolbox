%%% INPUT THE IMAGE FILE NAME:
if ~exist('fc','var')||~exist('cc','var')||~exist('kc','var')||~exist('alpha_c','var'),
    fprintf(1,'No intrinsic camera parameters available.\n');
    return;
end;

disp('Program that undistorts images:');
disp('The intrinsic camera parameters are assumed to be known (previously computed)');

if norm(kc)==0,
    fprintf(1,'The lens distortion seems to be zero!\n');
    return;
end;

fprintf(1,'\n');
quest = input('Do you want to undistort all the calibration images ([]=0) or a new image (1)? ');
if isempty(quest),
    quest = 0;
end;

if ~quest,      % rect all the images
    if n_ima == 0,
        fprintf(1,'No image data available\n');
        return;
    end;
    % READ IN IMAGE 1:
    ima_name = [calib_name strnum_cell{1} '.' format_image];
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
    [ny,nx,nc] = size(I);
    % Pre-compute the necessary indices and blending coefficients to enable quick rectification:
    [~,ind_new,ind_1,ind_2,ind_3,ind_4,a1,a2,a3,a4] = rect_index(zeros(ny,nx),eye(3),fc,cc,kc,alpha_c);
    
    for kk = 1:n_ima,
        ima_name = [calib_name strnum_cell{kk} '.' format_image];
        if  exist(ima_name,'file')~=2,
            fprintf(1,'Image %s not found!!!\n',ima_name);
        else
            fprintf(1,'Loading image %s...\n',ima_name);
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
            if nc>1,
                I = 0.299 * I(:,:,1) + 0.587 * I(:,:,2) + 0.114 * I(:,:,3);
            end;
            %    Irect = zeros(ny,nx);         % black background
            Irect = 255*ones(ny,nx);     % white background
            Irect(ind_new) = a1 .* I(ind_1) + a2 .* I(ind_2) + a3 .* I(ind_3) + a4 .* I(ind_4);
            
            ima_name2 = ['rect_' ima_name];
            fprintf(1,'Saving undistorted image as %s ...\n', ima_name2);
            if format_image(1) == 'p',
                if format_image(2) == 'p',
                    saveppm(ima_name2,uint8(round(Irect)));
                elseif format_image(2) == 'g',
                    savepgm(ima_name2,uint8(round(Irect)));
                else
                    imwrite(uint8(round(Irect)),gray(256),ima_name2,format_image);
                end;
            else
                if format_image(1) == 'r',
                    writeras(ima_name2,round(Irect),gray(256));
                else
                    imwrite(uint8(round(Irect)),gray(256),ima_name2,format_image);
                end;
            end;
        end;
    end;
    
else
    
    dir;
    fprintf(1,'\n');
    format_image2 = format_image;
    image_name = input('Image name (full name without extension): ','s');
    ima_name = [image_name '.' format_image2];
    if exist(ima_name,'file')~=2,
        fprintf(1,['\n' ima_name ' not found! Please re-input image format:\n']);
        format_image2 = '0';
        while format_image2 == '0',
            format_image2 =  input(['Image format: ( ''n''=''png'', ''b''=''bmp'', ''t''=''tif'', ' ...
                '''j''=''jpg'', ''g''=''jpeg'',\n''p''=''pgm'', ''m''=''ppm'', []=''r''=''ras'') '],'s');
            if isempty(format_image2),
                format_image2 = 'ras';
            end;
            if lower(format_image2(1)) == 'n',
                format_image2 = 'png';
            elseif lower(format_image2(1)) == 'b',
                format_image2 = 'bmp';
            elseif lower(format_image2(1)) == 't',
                format_image2 = 'tif';
            elseif lower(format_image2(1)) == 'j',
                format_image2 = 'jpg';
            elseif lower(format_image2(1)) == 'g',
                format_image2 = 'jpeg';
            elseif lower(format_image2(1)) == 'p',
                format_image2 = 'pgm';
            elseif lower(format_image2(1)) == 'm',
                format_image2 = 'ppm';
            elseif lower(format_image2(1)) == 'r',
                format_image2 = 'ras';
            else
                disp('Invalid image format');
                format_image2 = '0'; % Ask for format once again
            end;
        end;
        ima_name = [image_name '.' format_image2];
    end;
    % READ IN IMAGE:
    if format_image2(1) == 'p',
        if format_image2(2) == 'p',
            I = double(loadppm(ima_name));
        elseif format_image2(2) == 'g',
            I = double(loadpgm(ima_name));
        else
            I = double(imread(ima_name));
        end;
    else
        if format_image2(1) == 'r',
            I = readras(ima_name);
        else
            I = double(imread(ima_name));
        end;
    end;
    
    [m,n,junk] = size(I);
    if junk==3,
        I = 0.299 * I(:,:,1) + 0.5870 * I(:,:,2) + 0.114 * I(:,:,3);
    end;
    assert ((m==ny)&&(n==nx),'Image size do not match with the calibration data!');
    
    % SHOW THE ORIGINAL IMAGE:
    figure(2);
    image(I);
    colormap(gray(256));
    axis image;
    title('Original image (with distortion) - Stored in array I');
    drawnow;
    
    % UNDISTORT THE IMAGE:
    Irect = rectify_image(I,fc,cc,kc,alpha_c);
    figure(3);
    image(Irect);
    colormap(gray(256));
    axis image;
    title('Undistorted image - Stored in array Irect');
    drawnow;
    
    % SAVE THE IMAGE IN FILE:
    ima_name2 = ['rect_' ima_name];
    fprintf(1,['Saving undistorted image under ' ima_name2 '...']);
    if format_image2(1) == 'p',
        if format_image2(2) == 'p',
            saveppm(ima_name2,uint8(round(Irect)));
        elseif format_image2(2) == 'g',
            savepgm(ima_name2,uint8(round(Irect)));
        else
            imwrite(uint8(round(Irect)),gray(256),ima_name2,format_image2);
        end;
    else
        if format_image2(1) == 'r',
            writeras(ima_name2,round(Irect),gray(256));
        else
            imwrite(uint8(round(Irect)),gray(256),ima_name2,format_image2);
        end;
    end;
end;
fprintf(1,'Done!\n');