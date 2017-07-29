%  See also extract_foreground.
% This script extract foreground of image sequence using adaptive threshold.
% initialization of number of image frames
ndarea = 0.1;
n_frame = 0;
if exist('imgdir','var') && ~isempty(imgdir),
    if exist('imgbase','var') && ~isempty(imgbase) && exist('imgfmt','var') && ~isempty(imgfmt),
        %file name prefix number of original image
        imPreNum = [imgdir '/' imgbase];
        [n_frame, strnum_frame, frame_num] = check_image_sequence(imPreNum, imgfmt);
        if n_frame ==0,
            fprintf(1,'No images found! Please relocate!\n');
        end;
    end;
end;

while n_frame == 0,
    [imgbase,imgdir]= uigetfile({'*.jpg;*.bmp;*.png;*.gif;*.tif','Image Files (*.jpg,*.bmp,*.png,*.gif,*.tif)';...
        '*.jpg','JPEG Files(*.jpg)'; '*.bmp','BMP Files(*.bmp)'; '*.png','PNG Files(*.png)'; '*.gif','GIF Files(*.gif)';...
        '*.tif','TIFF Files(*.tif)'; '*.*','All Files(*.*)'},'Select the first frame image!');    % 选择图片路径
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
        [n_frame, strnum_frame, frame_num] = check_image_sequence(imPreNum, imgfmt);
    end;
end;

% skip step for long run check
if n_frame > 10,
    step = ceil(n_frame/10);       % sample 9 pics maximum
else
    step = 1;
end;
wt = 0.3;       % wait time for checking sample effect

flag = 0;    % switch to skip the main extraction process or not
% load variables of extraction setting if there is
save_name = [imgdir '/Extract_' imgbase '.mat'];
if exist(save_name,'file')==2,
    load(save_name);
    imPreNumf = [imgdir '/' fgprefix imgbase];
    if exist([imPreNumf strnum_frame{n_frame} '.png'],'file')==2,
        fprintf(1,'\nForeground binary images detected ...\n');
        flag = input('Skip the foregroud extracting process or not? ([]=yes, other=no) ','s');
        flag = isempty(flag);
    end;
else
    if ~exist('fgprefix','var') || isempty(fgprefix),
        fgprefix = input('Prefix of foreground images you want to save: ([]=BW_) ','s');
        if isempty(fgprefix),
            fgprefix = 'BW_';
        end;
    end;
    % file name prefix number of foreground image
    imPreNumf = [imgdir '/' fgprefix imgbase];
    interest_box = [];
end;

if flag,
    fprintf(1,'\nForeground extraction abort ...\n');
    % enter the upper folder of image directory
    cd([imgdir '/../']);
    return;
end;

% set active images
if ~exist('ind_active','var') || isempty(ind_active),
    ind_active = input(['Please input index number of active frames: ([]=1:' num2str(n_frame) ') ']);
    if isempty(ind_active),
        ind_active=1:n_frame;
    else
        ind_active = round(ind_active);
        assert(all((ind_active<=n_frame) & (ind_active>0)),'Unexpected image numbers!');
    end;
end;
N_active = length(ind_active);

% check image color channels (1 or 3)
frame_kk = imread([imPreNum strnum_frame{1} '.' imgfmt]);
[ny, nx, npage] = size(frame_kk);     % size of images
blob = ones(3);

flag = 0;    % switch to skip the foreground initialization process or not
if exist('interest_box','var') && ~isempty(interest_box),
    fprintf(1,'\nThe foreground interest box for refinement detected!\n');
    flag = input('Skip the foreground initialization process or not? ([]=yes, other=no) ','s');
    flag = isempty(flag);
end;
if ~flag,
    if ~exist('bgprefix','var') || isempty(bgprefix),
        bgprefix = input('Prefix of background images: ([]=back) ','s');
        if isempty(bgprefix),
            bgprefix = 'back';
        end;
    end;
    fprintf(1,'Make sure there are background images with prefix ''%s'' in the directory!\n', bgprefix);
    if ~exist('bgfmt','var') || isempty(bgfmt),
        bgfmt = input(['Format of background images: ([]=' imgfmt ') '],'s');
        if isempty(bgfmt),
            bgfmt  = imgfmt;
        end;
    end;
    bgname = [imgdir '/' bgprefix '_' imgbase '.' bgfmt];
    
    % load background
    if exist(bgname,'file')==2,
        fprintf(1,'\nLoading the image %s ...\n',[bgprefix '_' imgbase '.' bgfmt]);
        background = imread(bgname);
    else
        fprintf(1,'\nGenerating mean background image if there are multiple images with the specific prefix ...\n');
        % get the mean of several background images, and convert to gray
        background = uint8(imgrayMean(imgdir,bgprefix,bgfmt));
        if isempty(background),
            imgdir = [];
            fprintf(1,'\nERROR: You should prepare background images first!\n');
            return;
        end;
        imwrite(background,bgname);  % 保存灰度图像
    end;
    figure(2);
    image(background);
    colormap(gray(256));
    title('Image background');
    axis image;
    
    % initialization
    Not_save = 1;
    while Not_save,
         fprintf(1,'\nSet signal amplifier (>1) to enhance the difference between foreground and background ...\n');
        if exist('amplifier','var'),
            default_set = amplifier;
        else
            default_set = 5;
        end;
        amplifier = input(['amplifier = ([]=' num2str(default_set) ') ']);
        if isempty(amplifier),
            amplifier = default_set;
        end;
        
        fprintf(1,'\nSet the global threshold (0.0~1.0) to extract initial foreground ...\n');
        if exist('bw_thresh','var'),
            default_set = bw_thresh;
        else
            default_set = 0.1;
        end;
        bw_thresh = input(['bw_thresh = ([]=' num2str(default_set) ') ']);
        if isempty(bw_thresh),
            bw_thresh = default_set;
        end;
        
        % check for the long run, sample at large step (check 10 pictures maximum)
        for kk = ind_active(1 : step : N_active),   % no need to display frame1 again
            frame_kk = imread([imPreNum strnum_frame{kk} '.' imgfmt]);
            if npage==3,
                frame_kk = 0.299 * frame_kk(:,:,1) + 0.587 * frame_kk(:,:,2) + 0.114 * frame_kk(:,:,3);
            end;
            % 背景减去帧序列，得到白色的昆虫和黑色的环境 (imsubtract)
            fore_kk = (background - frame_kk)*amplifier;
            fore_kk = imclose(im2bw(fore_kk, bw_thresh),blob);  % Thresholds to get foreground
            % show the effect of sampled frames
            figure(2);
            image(fore_kk);
            colormap(gray(2));      % colormap for binary image
            title(['Check initial foreground of frame: ' num2str(kk)]);
            axis image;
            pause(wt);           % wait for wt second
        end;
        fprintf(1,'\nCheck the extracted foreground ...\n');
        Not_save = input('Need to reset the global threshold? ([]=yes, other=no) ','s');
        Not_save = isempty(Not_save);
    end;
    
    %%% delete regions of error forground
    del_region = false(0); nn = 0;
    re_edit = input('Need to delete some trash region or not? ([]=yes, other=no) ','s');
    re_edit = isempty(re_edit);
    if re_edit,
        % generate meshgrid for inpolygon
        [xxg, yyg] = meshgrid(1:nx, 1:ny);
    end;
    while re_edit,
        Not_save = 1;
        while Not_save,
            frame_kk = imread([imPreNum strnum_frame{ind_active(1)} '.' imgfmt]);
             if npage==3,
                frame_kk = 0.299 * frame_kk(:,:,1) + 0.587 * frame_kk(:,:,2) + 0.114 * frame_kk(:,:,3);
            end;
            fore_kk = (background - frame_kk)*amplifier;
            fore_kk = imclose(im2bw(fore_kk, bw_thresh),blob);
            for i = 1:nn,
                fore_kk = fore_kk & ~del_region(:,:,i);
            end;
            fprintf(1,['\nPlease click polygon vertices with left mouse button first\n' ...
                'and click the last point with right mouse button ...\n']);
            figure(2);
            image(fore_kk);
            colormap(gray(2));      % colormap for binary image
            title('LMB to choose vertices - RMB to finish clicking:');
            axis image;
            hold on;
            % choose vertices for polygon
            [xp,yp] = ginput(1);
            xvert = xp;
            yvert = yp;
            plot(xp,yp,'go');
            lmb = 1;      % ginput "left mouse button": 1=LMB, 2=MMB, 3=RMB
            while lmb==1 || length(xvert)<3,
                [xp,yp,lmb] = ginput(1);
                xvert = [xvert, xp];
                yvert = [yvert, yp];
                plot(xp,yp,'go');
                % do not cover old plot
                plot(xvert(end-1:end),yvert(end-1:end),'b-','linewidth',2);
            end;
            plot([xp,xvert(1)], [yp,yvert(1)],'b-','linewidth',2);
            hold off;
            % check the polygon in samples
            for kk = ind_active(step+1 : step : N_active),    % no need to display frame1 again
                frame_kk = imread([imPreNum strnum_frame{kk} '.' imgfmt]);
                if npage==3,
                    frame_kk = 0.299 * frame_kk(:,:,1) + 0.587 * frame_kk(:,:,2) + 0.114 * frame_kk(:,:,3);
                end;
                fore_kk = (background - frame_kk)*amplifier;
                fore_kk = imclose(im2bw(fore_kk, bw_thresh),blob);
                for i = 1:nn,
                    fore_kk = fore_kk & ~del_region(:,:,i);
                end;
                % show the effect of sampled frames
                figure(2);
                image(fore_kk);
                colormap(gray(2));      % colormap for binary image
                title('Check the polygon in sample images:');
                axis image;
                hold on;
                plot([xvert, xvert(1)], [yvert, yvert(1)],'b-','linewidth',2);
                hold off;
                pause(wt);           % wait for wt second
            end;
            Not_save = input('Need to readjust the polygon region for deleting? ([]=yes, other=no) ','s');
            Not_save = isempty(Not_save);
        end;
        % binary image of polygon in region
        del_region = cat(3,del_region, inpolygon(xxg, yyg, xvert, yvert));
        nn = nn+1;
        % display results after amendment
        for kk = ind_active(1 : step : N_active),
            frame_kk = imread([imPreNum strnum_frame{kk} '.' imgfmt]);
             if npage==3,
                frame_kk = 0.299 * frame_kk(:,:,1) + 0.587 * frame_kk(:,:,2) + 0.114 * frame_kk(:,:,3);
            end;
            fore_kk = (background - frame_kk)*amplifier;
            fore_kk = imclose(im2bw(fore_kk, bw_thresh),blob);
            for i = 1:nn,
                fore_kk = fore_kk & ~del_region(:,:,i);
            end;
            % show the effect of sampled frames
            figure(2);
            image(fore_kk);
            colormap(gray(2));      % colormap for binary image
            title('Check foreground images after amendment:');
            axis image;
            pause(wt);           % wait for wt second
        end;
        re_edit = input('Still need amendment or not? ([]=yes, other=no) ','s');
        re_edit = isempty(re_edit);
    end;
    
    % compute interest box for every image: subscipts of left-up corner and right-bottom corner.
    interest_box = cell(1,n_frame);
    for kk = ind_active,
        frame_kk = imread([imPreNum strnum_frame{kk} '.' imgfmt]);
        if npage==3,
            frame_kk = 0.299 * frame_kk(:,:,1) + 0.587 * frame_kk(:,:,2) + 0.114 * frame_kk(:,:,3);
        end;
        fore_kk = (background - frame_kk)*amplifier;
        fore_kk = imclose(im2bw(fore_kk, bw_thresh),blob);
        for i = 1:nn,
            fore_kk = fore_kk & ~del_region(:,:,i);
        end;
        geom_kk = regionprops(fore_kk);
        geom_kk = geom_kk([geom_kk.Area]>sum([geom_kk.Area])*ndarea);
        n = length(geom_kk);
        temp = cat(1,geom_kk.BoundingBox)';
        temp(3:4,:) = min(temp([2,1],:)+temp([4,3],:)+9.5, [ny; nx]*ones(1,n));
        temp(1:2,:) = max(temp([2,1],:)-9.5,ones(2,n));
        interest_box{kk} = temp;
    end;
end;

background = imread([imgdir '/' bgprefix '_' imgbase '.' bgfmt]);
nn = size(del_region,3);
%%% refine foreground by means of Canny edge and adaptive threshold
Not_save = 1;
while Not_save,
    fprintf(1,'\nSet parameters for the anisotropic diffusion function to enhance image:\n');
    fprintf(1,'\nSet the gradient modulus kappa (>0) for the anisotropic diffusion function ...\n');
    if exist('kappa','var'),
        default_set = kappa;
    else
        default_set = 3;
    end;
    kappa = input(['kappa = ([]=' num2str(default_set) ') ']);
    if isempty(kappa),
        kappa = default_set;
    end;
    
    fprintf(1,'\nSet the iteration number (1~30) for the anisotropic diffusion function ...\n');
    if exist('num_iter','var'),
        default_set = num_iter;
    else
        default_set = 10;
    end;
    num_iter = input(['num_iter = ([]=' num2str(default_set) ') ']);
    if isempty(num_iter),
        num_iter = default_set;
    end;
    
    fprintf(1,'\nSet the switch of conduction coefficient for the anisotropic diffusion function ...\n');
    if exist('privil_sw','var'),
        default_set = privil_sw;
    else
        default_set = 1;
    end;
    fprintf(1,['\n1: privileges wide regions over smaller ones.' ...
        '\n0: privileges high-contrast edges over low-contrast ones.\n']);
    privil_sw = input(['Do you want to privilege wide regions over smaller ones? ([]=' num2str(default_set) ') ']);
    if isempty(privil_sw),
        privil_sw = default_set;
    else
        privil_sw = ~~privil_sw;
    end;
    
    delta_t = 1/7;
    % check for the long run, sample at large step (check 10 pictures maximum)
    for kk = ind_active(1 : step : N_active),
        frame_kk = imread([imPreNum strnum_frame{kk} '.' imgfmt]);
        if npage==3,
            frame_kk = 0.299 * frame_kk(:,:,1) + 0.587 * frame_kk(:,:,2) + 0.114 * frame_kk(:,:,3);
        end;
        frame_kk = background - frame_kk;
        for i = 1:nn,
            frame_kk(del_region(:,:,i)) = 0;
        end;
        ibox = interest_box{kk};
        n = size(ibox,2);
        figure(2);
        for i=1:n,
            temp = frame_kk(ibox(1,i) : ibox(3,i), ibox(2,i) : ibox(4,i));
            temp =  anisodiff2D(temp, kappa, num_iter, privil_sw, delta_t)*amplifier;
            subplot(1,n,i);
            image(temp);
            axis image;
        end;
        colormap(gray(256));      % colormap for grayscale image
        % title(['Enhanced foreground of frame: ', num2str(kk)]);
        pause(wt);           % wait for wt second
    end;
    fprintf(1,'\nCheck the enhanced foreground ...\n');
    Not_save = input('Need to reset the parameters? ([]=yes, other=no) ','s');
    Not_save = isempty(Not_save);
end;

Not_save = 1;
while Not_save,
    fprintf(1,'\nSet signal amplifier (>1) again to extract foreground from the enhanced image ...\n');
    default_set = amplifier;
    amplifier = input(['amplifier = ([]=' num2str(default_set) ') ']);
    if isempty(amplifier),
        amplifier = default_set;
    end;
    
    fprintf(1,'\nSet global threshold (0.0~1.0) again to extract foreground from enhanced image:\n');
    default_set = bw_thresh;
    bw_thresh = input(['bw_thresh = ([]=' num2str(default_set) ') ']);
    if isempty(bw_thresh),
        bw_thresh = default_set;
    end;

    fprintf(1,'\nSet the switch to patch holes in foreground (Be cautious with 1=yes) ... \n');
    if exist('patch_sw','var'),
        default_set = patch_sw;
    else
        default_set = 0;
    end;
    patch_sw = input(['Do you want to patch holes in foreground? ([]=' num2str(default_set) ') ']);
    if isempty(patch_sw),
        patch_sw = default_set;
    else
        patch_sw = ~~patch_sw;
    end;

    % check for the long run, sample at large step (check 10 pictures maximum)
    for kk = ind_active(1 : step : N_active),
        frame_kk = imread([imPreNum strnum_frame{kk} '.' imgfmt]);
        if npage==3,
            frame_kk = 0.299 * frame_kk(:,:,1) + 0.587 * frame_kk(:,:,2) + 0.114 * frame_kk(:,:,3);
        end;
        frame_kk = background - frame_kk;
        for i = 1:nn,
            frame_kk(del_region(:,:,i)) = 0;
        end;
        ibox = interest_box{kk};
        n = size(ibox,2);
        figure(2);
        for i=1:n,
            temp = frame_kk(ibox(1,i) : ibox(3,i), ibox(2,i) : ibox(4,i));
            temp =  anisodiff2D(temp, kappa, num_iter, privil_sw, delta_t)*amplifier;
            fore_kk = imclose(im2bw(uint8(temp), bw_thresh),blob);
            if patch_sw,
                fore_kk = imfill(fore_kk, 'holes');
            end;
            subplot(1,n,i);
            image(fore_kk);
            axis image;
        end;
        colormap(gray(2));      % colormap for binary image
        % title(['Foreground edge of frame: ', num2str(kk)]);
        pause(wt);           % wait for wt second
    end;
    fprintf(1,'\nCheck the extracted foreground ...\n');
    Not_save = input('Need to reset the parameters? ([]=yes, other=no) ','s');
    Not_save = isempty(Not_save);
end;

Not_save = 1;
while Not_save,
    fprintf(1,'\nSet parameters for adaptive threshold function to refine foreground:\n');
    fprintf(1,'\nSet two gaussian filter (0~25) with different size for the adaptive threshold function ...\n');
    if exist('gauss_sig','var'),
        default_set = gauss_sig;
    else
        default_set = [10, 3];
    end;
    gauss_sig = input(['gauss_sig = ([]=[' num2str(default_set) ']) ']);
    if isempty(gauss_sig),
        gauss_sig = default_set;
    end;
    if length(gauss_sig)==1,
        gauss_sig(2) = 0.3*gauss_sig;
    end;
    
    fprintf(1,'\nSet the fixed threshold (-20~20) for the adaptive threshold function ...\n');
    if exist('fixed_thresh','var'),
        default_set = fixed_thresh;
    else
        default_set = -1;
    end;
    fixed_thresh = input(['fixed_thresh = ([]=' num2str(default_set) ') ']);
    if isempty(fixed_thresh),
        fixed_thresh = default_set;
    end;

    % check for the long run, sample at large step (check 10 pictures maximum)
    for kk = ind_active(1 : step : N_active),
        frame_kk = imread([imPreNum strnum_frame{kk} '.' imgfmt]);
        if npage==3,
            frame_kk = 0.299 * frame_kk(:,:,1) + 0.587 * frame_kk(:,:,2) + 0.114 * frame_kk(:,:,3);
        end;
        frame_kk = background - frame_kk;
        for i = 1:nn,
            frame_kk(del_region(:,:,i)) = 0;
        end;
        ibox = interest_box{kk};
        n = size(ibox,2);
        figure(2);
        for i=1:n,
            temp = frame_kk(ibox(1,i) : ibox(3,i), ibox(2,i) : ibox(4,i));
            temp =  anisodiff2D(temp, kappa, num_iter, privil_sw, delta_t)*amplifier;
            BW = imclose(im2bw(uint8(temp), bw_thresh),blob);
            if patch_sw,
                BW = imfill(BW, 'holes');
            end;
            for j=1:2,
                fore_kk = BW & adaptivethresh(temp, gauss_sig(j), fixed_thresh, 'gaussian','fixed');
                subplot(2,n, i+(j-1)*n);
                image(fore_kk);
                axis image;
            end;
        end;
        colormap(gray(2));      % colormap for binary image
        % title(['Foreground edge of frame: ', num2str(kk)]);
        pause(wt);           % wait for wt second
    end;
    fprintf(1,'\nCheck the extracted foreground ...\n');
    Not_save = input('Need to reset the parameters? ([]=yes, other=no) ','s');
    Not_save = isempty(Not_save);
end;

close(2);
Not_save = 1;
while Not_save,
    fprintf(1,'\nSet the switch to fill foreground edges ...\n');
    if exist('filledge_sw','var'),
        default_set = filledge_sw;
    else
        default_set = 0;
    end;
    filledge_sw = input(['Extract and fill foreground edges or not? ([]=' num2str(default_set) ') ']);
    if isempty(filledge_sw),
        filledge_sw = default_set;
    else
        filledge_sw = ~~filledge_sw;
    end;
    
    if filledge_sw,
        fprintf(1,'\nSet parameters to fill gaps of foreground edges:\n');
        fprintf(1,'\nSet the global window size (>=4) for closeEdge function ...\n');
        if exist('gap_radius','var'),
            default_set = gap_radius;
        else
            default_set = 10;
        end;
        gap_radius = input(['gap_radius = ([]=' num2str(default_set) ') ']);
        if isempty(gap_radius),
            gap_radius = default_set;
        end;
        
        fprintf(1,'\nSet the switch to delete short edges for closeEdge function ...\n');
        if exist('dedge_sw','var'),
            default_set = dedge_sw;
        else
            default_set = 1;
        end;
        fprintf(1,['\n1: Delete short edges and fill gaps all at once.' ...
            '\n0: Keep short edges and fill gaps iteratively.\n']);
        dedge_sw = input(['Do you want to delete short edges? ([]=' num2str(default_set) ') ']);
        if isempty(dedge_sw),
            dedge_sw = default_set;
        else
            dedge_sw = ~~dedge_sw;
        end;
    end;
    
    fprintf(1,'\nSet the switch to delete small blobs (as noise) in foreground at last ...\n');
    if exist('dblob_sw','var'),
        default_set = dblob_sw;
    else
        default_set = 1;
    end;
    fprintf(1,'\nBlobs of whose area less than %g of the biggest 4-connected domain will be deleted!\n', ndarea);
    dblob_sw = input(['Delete small blobs or not? ([]=' num2str(default_set) ') ']);
    if isempty(dblob_sw),
        dblob_sw = default_set;
    else
        dblob_sw = ~~dblob_sw;
    end;
    
    % check for the long run, sample at large step (check 10 pictures maximum)
    for kk = ind_active(1 : step : N_active),
        frame_kk = imread([imPreNum strnum_frame{kk} '.' imgfmt]);
        if npage==3,
            frame_kk = 0.299 * frame_kk(:,:,1) + 0.587 * frame_kk(:,:,2) + 0.114 * frame_kk(:,:,3);
        end;
        frame_kk = background - frame_kk;
        for i = 1:nn,
            frame_kk(del_region(:,:,i)) = 0;
        end;
        fore_kk = false(ny,nx);
        ibox = interest_box{kk};
        n = size(ibox,2);
        for i=1:n,
            temp = frame_kk(ibox(1,i) : ibox(3,i), ibox(2,i) : ibox(4,i));
            temp =  anisodiff2D(temp, kappa, num_iter, privil_sw, delta_t)*amplifier;
            BW =  imclose(im2bw(uint8(temp), bw_thresh),blob);
             if patch_sw,
                BW = imfill(BW, 'holes');
            end;
            BW1 = BW & (adaptivethresh(temp, gauss_sig(1), fixed_thresh, 'gaussian','fixed') ...
                | adaptivethresh(temp, gauss_sig(2), fixed_thresh, 'gaussian','fixed'));
            if filledge_sw,
                BW2 = closeEdge(canny_edge(temp), gap_radius, dedge_sw);
                % keep holes in foreground
                BW1 = BW & imfill(BW1 | BW2, 'holes');
            end;
            if dblob_sw,        % erase small blobs in the foreground
                temp = bwconncomp(BW1,4);
                numPixels = cellfun(@numel,temp.PixelIdxList);
                idx = numPixels<max(numPixels)*ndarea;
                BW1(cell2mat({temp.PixelIdxList{idx}}')) = 0;
            end;
            fore_kk(ibox(1,i) : ibox(3,i), ibox(2,i) : ibox(4,i)) = BW1;
        end;
        % show the effect of sampled frames
        figure(2);
        image(fore_kk);
        colormap(gray(2));      % colormap for binary image
        title(['Refined foreground of frame: ', num2str(kk)]);
        axis image;
        pause(wt);           % wait for wt second
    end;
    fprintf(1,'\nCheck the refined foreground ...\n');
    Not_save = input('Need to reset the parameters? ([]=yes, other=no) ','s');
    Not_save = isempty(Not_save);
end;

% save refined foreground images
fprintf(1,'\nSaving refined foreground for every frame, please wait ...\n');
for kk = ind_active,
    frame_kk = imread([imPreNum strnum_frame{kk} '.' imgfmt]);
    if npage==3,
        frame_kk = 0.299 * frame_kk(:,:,1) + 0.587 * frame_kk(:,:,2) + 0.114 * frame_kk(:,:,3);
    end;
    frame_kk = background - frame_kk;
    for i = 1:nn,
        frame_kk(del_region(:,:,i)) = 0;
    end;
    fore_kk = false(ny,nx);
    ibox = interest_box{kk};
    n = size(ibox,2);
    for i=1:n,
        temp = frame_kk(ibox(1,i) : ibox(3,i), ibox(2,i) : ibox(4,i));
        temp =  anisodiff2D(temp, kappa, num_iter, privil_sw, delta_t)*amplifier;
        BW =  imclose(im2bw(uint8(temp), bw_thresh),blob);
        if patch_sw,
            BW = imfill(BW, 'holes');
        end;
        BW1 = BW & (adaptivethresh(temp, gauss_sig(1), fixed_thresh, 'gaussian','fixed') ...
            | adaptivethresh(temp, gauss_sig(2), fixed_thresh, 'gaussian','fixed'));
        if filledge_sw,
            BW2 = closeEdge(canny_edge(temp), gap_radius, dedge_sw);
            % keep holes in foreground
            BW1 = BW & imfill(BW1 | BW2, 'holes');
        end;
        if dblob_sw,        % erase small blobs in the foreground
            temp = bwconncomp(BW1,4);
            numPixels = cellfun(@numel,temp.PixelIdxList);
            idx = numPixels<max(numPixels)*ndarea;
            BW1(cell2mat({temp.PixelIdxList{idx}}')) = 0;
        end;
        fore_kk(ibox(1,i) : ibox(3,i), ibox(2,i) : ibox(4,i)) = BW1;
    end;
    imwrite(fore_kk,[imPreNumf strnum_frame{kk} '.png'],'png');   % background设为0
end;
fprintf(1,'\nImage initial foreground all saved!\n');

% save the variable settings for extration
string_save = ['save ' save_name ' bgprefix bgfmt fgprefix amplifier bw_thresh patch_sw del_region ' ...
    'interest_box kappa num_iter privil_sw gauss_sig fixed_thresh filledge_sw dblob_sw'];
if filledge_sw,
    string_save = [string_save ' gap_radius dedge_sw'];
end;
eval(string_save);

% enter the upper folder of image directory
cd([imgdir '/../']);