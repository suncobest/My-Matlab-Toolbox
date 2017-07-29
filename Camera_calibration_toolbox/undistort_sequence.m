%%% INPUT THE IMAGE FILE NAME:

if ~exist('drawIM','var'),
    drawIM = input('Display the original and undistorted images or not? ([]=no, other=yes) ','s');
    drawIM = ~isempty(drawIM);
end;

if ~exist('fc','var')||~exist('cc','var')||~exist('kc','var')||~exist('alpha_c','var'),
   fprintf(1,'No intrinsic camera parameters available. Maybe, need to load Calib_Results*.mat\n');
   return;
end;

disp('Program that undistorts a whole sequence of images.');
disp('The intrinsic camera parameters are assumed to be known (previously computed)!');

dir;
fprintf(1,'\n');
seq_name = input('Basename of sequence images (without number nor suffix): ','s');

format_imseq = '0';
while format_imseq == '0',
    format_imseq =  input(['Image format: ( ''n''=''png'', ''b''=''bmp'', ''t''=''tif'', ' ...
        '''j''=''jpg'', ''g''=''jpeg'',\n''p''=''pgm'', ''m''=''ppm'', []=''r''=''ras'') '],'s');
    if isempty(format_imseq),
        format_imseq = 'ras';
    end;
    if lower(format_imseq(1)) == 'n',
        format_imseq = 'png';
    elseif lower(format_imseq(1)) == 'b',
        format_imseq = 'bmp';
    elseif lower(format_imseq(1)) == 't',
        format_imseq = 'tif';
    elseif lower(format_imseq(1)) == 'j',
        format_imseq = 'jpg';
    elseif lower(format_imseq(1)) == 'g',
        format_imseq = 'jpeg';
    elseif lower(format_imseq(1)) == 'p',
        format_imseq = 'pgm';
    elseif lower(format_imseq(1)) == 'm',
        format_imseq = 'ppm';
    elseif lower(format_imseq(1)) == 'r',
        format_imseq = 'ras';
    else
        disp('Invalid image format');
        format_imseq = '0'; % Ask for format once again
    end;
end;

ima_sequence = dir([seq_name '*.' format_imseq]);
if isempty(ima_sequence),
    fprintf(1,'No image found\n');
    return;
end;

ima_name = ima_sequence(1).name;
if format_imseq(1) == 'p',
    if format_imseq(2) == 'p',
        I = double(loadppm(ima_name));
    elseif format_imseq(2) == 'g',
        I = double(loadpgm(ima_name));
    else
        I = double(imread(ima_name));
    end;
else
    if format_imseq(1) == 'r',
        I = readras(ima_name);
    else
        I = double(imread(ima_name));
    end;
end;
[ny,nx,nc] = size(I);

% Pre-compute the necessary indices and blending coefficients to enable quick rectification: 
[~,ind_new,ind_1,ind_2,ind_3,ind_4,a1,a2,a3,a4] = rect_index(zeros(ny,nx),eye(3),fc,cc,kc,alpha_c);
n_seq = length(ima_sequence);

for kk = 1:n_seq,
    ima_name = ima_sequence(kk).name;
    fprintf(1,'Loading original image %s...',ima_name);
    
    %%% READ IN IMAGE:
    if format_imseq(1) == 'p',
        if format_imseq(2) == 'p',
            I = double(loadppm(ima_name));
        elseif format_imseq(2) == 'g',
            I = double(loadpgm(ima_name));
        else
            I = double(imread(ima_name));
        end;
    else
        if format_imseq(1) == 'r',
            I = readras(ima_name);
        else
            I = double(imread(ima_name));
        end;
    end;
    
%     Irect = zeros(ny,nx,nc);     % black background
    Irect = 255*ones(ny,nx,nc);     % white background
    for ii = 1:nc,
        Iii = I(:,:,ii);
        I2ii = Irect(:,:,ii);
        I2ii(ind_new) = a1 .* Iii(ind_1) + a2 .* Iii(ind_2) + a3 .* Iii(ind_3) + a4 .* Iii(ind_4);
        Irect(:,:,ii) = I2ii;
    end;
    Irect = uint8(Irect);
    if drawIM,
        figure(2);
        image(uint8(I));
        if nc==1,
            colormap(gray(256));
        end;
        drawnow;
        figure(3);
        image(Irect);
        if nc==1,
            colormap(gray(256));
        end;
        drawnow;
    end;
    ima_name2 = ['undist_' ima_name];
    fprintf(1,'Saving undistorted image under %s...\n',ima_name2);
    if format_imseq(1) == 'p',
        if format_imseq(2) == 'p',
            saveppm(ima_name2,Irect));
        elseif ormat_imseq(2) == 'g',
            savepgm(ima_name2,Irect));
        else
            imwrite(Irect,ima_name2,format_imseq);
        end;
    else
        if format_imseq(1) == 'r',
            writeras(ima_name2,Irect,gray(256));
        else
            imwrite(Irect,ima_name2,format_imseq);
        end;
    end;
end;
