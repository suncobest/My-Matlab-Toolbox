% find the codedir of the toolbox
codeDir = mfilename('fullpath');
flag = find(codeDir=='\' | codeDir=='/',1,'last');
if ~isempty(flag),
    codeDir = codeDir(1 : flag(end)-1);
end;
if exist([codeDir '/load_imgdir.m'],'file'),
    load_imgdir;
    fprintf(1,'\nModify file ''load_imgdir.m'' in code directory if you want to redirect image path!\n');
else
    fprintf(1,'\n''You should run scripts ''load_cameras'' and ''foreground_extraction'' first.\n');
    return;
end;

if ~exist('n_frame','var'),
    save_name = [imgdir '/Extraction_result.mat'];
    if exist(save_name,'file')==2,
        load(save_name);
    else
        fprintf(1,'\nYou should click button ''Foreground extraction'' first ...\n');
        return;
    end;
end;

if ~exist('Tcc','var'),
    save_name = [imgdir '/Camera_list.mat'];
    if exist(save_name,'file')==2,
        load(save_name);
    else
        fprintf(1,['\nYou should copy ''Camera_list.mat'' to the foreground ' ...
            'images'' directory!\nLocate the file ...\n']);
        camfolder = uigetdir;
        save_name = [camfolder '/Camera_list.mat'];
        if exist(save_name,'file')==2,
            load(save_name);
            copyfile(save_name, [imgdir '/Camera_list.mat']);
        else
            fprintf(1,'\n''Camera_list.mat'' not found!\n');
            return;
        end;
    end;
end;

save_name = [imgdir '/visualhull_environment.mat'];
if exist(save_name,'file')==2,
    fprintf(1,'\n''visualhull_environment.mat'' detected! You can load the file or start a new mission.\n');
    flag = input('Load ''visualhull_environment.mat'' or not: ([]=yes, other=no) ','s');
    if isempty(flag),
        fprintf(1,'\nLoading file ''visualhull_environment.mat'' ...\n');
        load(save_name);
        if level==1,
            fprintf(1,'Visual hull of level 1 detected, You can refine it now ...\n');
        else
            fprintf(1,'Visual hull of level %d detected, please check if reconstruction is complete ...\n', level);
        end;
        return;
    end;
end;

% check if foreground images are there
if ~exist('active_imgviews','var') || size(active_imgviews,2)~=n_frame,
    active_imgviews = true(n_cam, n_frame);
end;
if calib_mode,
    imPreNumf = [imgdir '/' fgprefix imgbase];
    ind_active = find(sum(active_imgviews,1)>=2);
    for kk=ind_active,
        if exist([imPreNumf strnum_frame{kk} '.png'],'file')~=2,
            fprintf(1,'\nForeground image of frame %d is missing, setting all views inactive!\n',kk);
            active_imgviews(:,kk) = 0;
        end;
    end;
else
    for pp = 1:n_cam,
        imgdir = imgdir_cell{pp};
        imgbase = imgbase_cell{pp};
        strnum_frame = strnum_frameNcam{pp};
        imPreNumf = [imgdir '/' fgprefix imgbase];
        ind_active = find(active_imgviews(pp,:));
        for kk = ind_active,
            if exist([imPreNumf strnum_frame{kk} '.png'],'file')~=2,
                fprintf(1,'\nForeground image of  (camera %d, frame %d) is missing, setting the image inactive!\n',pp,kk);
                active_imgviews(pp,kk) = 0;
            end;
        end;
    end;
end;
active_images = sum(active_imgviews,1)>=2;
ind_active = find(active_images);
if isempty(ind_active),
    fprintf(1,'\nERROR: It takes at least two camera views to caculate 3d! No frame is available!\n');
    return;
end;


%% initialization of 3D reconstruction

fprintf(1,['\nThis script will initialize visual hull computation.\n' ...
    'Make sure foreground image sequence are now distortion free!\n']);
if ~exist('character_len','var'),
    fprintf(1,'\nPlease input the characteristic length of 3D object you want to reconstruct ...\n');
    character_len = input('character_len = ([]=10) ');
    if isempty(character_len),
        character_len = 10;    % bee's wing length (unit: mm)
    end;
end;
% 'nxyz * character_len' equals the size of 3D bounding box; nxyz =[nla;nlb;nlc];
if ~exist('nxyz','var') || length(nxyz)~=3,
    fprintf(1,['\nPlease input the size of bounding box for reconstruction ' ...
        'in unit of characteristic length ...\n']);
    flag = 1;
    while flag,
        nxyz = input('nxyz = ([]=[3;3;3]) ');
        if isempty(nxyz),
            nxyz = [3;3;3];
            flag=0;
        elseif length(nxyz)==1,
            nxyz = nxyz(ones(3,1));
            flag=0;
        else
            nxyz = nxyz(:);
            flag = length(nxyz)~=3;
            if flag,
                fprintf(1,'Unexpected input!');
            end;
        end;
    end;
end;
boundsize = character_len*nxyz;

if ~exist('ndarea','var'),
    fprintf(1,['\nIf the foreground area of a view is small than a fraction of the 1st frame,\n' ...
        'then the camera view will be set inactive! Please input the fraction number...\n']);
    ndarea = input('Threshold fraction number = ([]=0.5) ');
    if isempty(ndarea),
        ndarea = 0.5;
    else
        assert(ndarea>0 && ndarea<1, 'Unexpected input!');
    end;
end;

if ~exist('tol_factor','var'),
    fprintf(1,['\nIf the next center position is too far from current one,\n' ...
        'then the camera view will be set inactive!\nPlease input value of ' ...
        'distance tolerance in unit of pixel bounding box ...\n']);
    tol_factor = input('Factor of distance tolerance = ([]=1 bounding box) ');
    if isempty(tol_factor) || tol_factor <=0,
        tol_factor = 1;
    end;
end;

saveIM = input('Draw and save 2D, 3D images or not? ([]=no, other=yes) ','s');
saveIM = ~isempty(saveIM);

n_view = n_frame*n_cam;
if ~exist('bounding_mat','var'),
    bounding_mat = NaN(4,n_view);   % roi: region of interest
    center2d_mat = NaN(2,n_view);
    area_mat = NaN(1,n_view);
    foreground_cell = cell(1,n_view);   % all foreground
end;

% image height in inch. 1pound = 1/72inch, 2inch can hold 12 lines of 12 pound words
height = 3;
if calib_mode,
    fprintf(1,'\nCompute image geometric feature under one-camera-with-split-views configuration!\n');
    kk = ind_active(1);
    active_view = active_imgviews(:,kk);
    active_imgviews = active_imgviews & active_view(:,ones(1,n_frame));
    ind = find(~active_view);
    for pp = ind,
      fprintf(1,'\nWarning: Camera view %d is inactive\n.',pp);
    end;
    frame_kk = imread([imPreNumf strnum_frame{kk} '.png']);
    % initialization
    [y,x] = find(frame_kk);
    xx = [x,y];
    nc = sum(active_view);
    centroid_last = NaN(n_cam,2);
    distance_toler = NaN(1,n_cam);
    threshold_area = distance_toler;
    while true,
        fprintf(1,'\nDistinguish camera views from each other:\n');
        figure(2); hold off;
        image(frame_kk);
        colormap(gray(2));      % colormap for binary image
        axis image;
        hold on;
        for pp = 1:n_cam,
            if active_view(pp),
                title(['Please click the foreground center of view ' num2str(pp) ':']);
                x = ginput(1);
                plot(x(1),x(2),'r+');
                centroid_last(pp,:) = x;
            end;
        end;
        [idx,center] = kmeans(xx,nc,'start',centroid_last(active_view,:));   % plant seeds
        for pp = 1:n_cam,
            if active_view(pp),
                % find the closest center to the clicked position
                [~, id] = min(sum((center-repmat(centroid_last(pp,:),nc,1)).^2, 2), [], 1);
                centroid_kk = center(id,:);
                ind = idx==id;
                area_kk = sum(ind);
                x = xx(ind,:)';
                temp = [min(x,[],2); max(x,[],2)];
                boundbox_kk = [temp(1:2)-0.5; temp(3:4)-temp(1:2)+1];
                distance_toler(pp) = max(boundbox_kk(3:4))*tol_factor;      % the longer edge of boundingbox
                threshold_area(pp) = area_kk*ndarea;
                kth = (kk-1)*n_cam+pp;
                area_mat(kth) = area_kk;
                bounding_mat(:,kth) = boundbox_kk;
                centroid_last(pp,:) = centroid_kk;
                center2d_mat(:,kth) = centroid_kk';
                foreground_cell{kth} = frame_kk(temp(2) : temp(4), temp(1) : temp(3));
                % plot camera number
                text(centroid_kk(1),centroid_kk(2),['\it\fontname{Arial}\fontsize{9}\color{green}' ...
                    num2str(pp)],'HorizontalAlignment','center');
                rectangle('Position',boundbox_kk,'edgecolor','b','linewidth',1);
            end;
        end;
        title('Check number of all camera views:');
        flag = input('Need to change camera view number or not? ([]=no, other=yes) ','s');
        if isempty(flag),
            break;
        end;
    end;
    set(gca,'position',[0 0 1 1]);
    resolution = round(ny/height);   % dpi
    set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
    drawnow;
    if saveIM,
        save_name = [imgdir '/loc_' imgbase strnum_frame{kk} '.jpg'];
        print(2,'-djpeg',['-r' num2str(resolution)],save_name);
    end;
    % load geometric feature of all image frames
    for kk = ind_active(2:end),
        frame_kk = imread([imPreNumf strnum_frame{kk} '.png']);
        active_view = active_imgviews(:,kk);
        % geometric feature
        [y,x] = find(frame_kk);
        xx = [x,y];
        [idx,center] = kmeans(xx,nc,'start',centroid_last(active_view,:));    % plant seeds
        figure(2); hold off;
        image(frame_kk);
        colormap(gray(2));      % colormap for binary image
        axis image;
        hold on;
        for pp = 1:n_cam,
            if active_view(pp),
                % find the closest center to center of last frame
                [y, id] = min(sum((center-repmat(centroid_last(pp,:),nc,1)).^2, 2), [], 1);
                centroid_kk = center(id,:);
                ind = idx==id;
                area_kk = sum(ind);
                x = xx(ind,:)';
                temp = [min(x,[],2); max(x,[],2)];
                boundbox_kk = [temp(1:2)-0.5; temp(3:4)-temp(1:2)+1];
                if sqrt(y)>distance_toler(pp) || area_kk<threshold_area(pp),
                  fprintf(1,['\nWarning: The foreground of (view %d, frame %d) is abnormal! '...
                               'Lost tracking of view %d...\n'],pp,kk,pp);
                    active_view(pp) = 0;
                    nc = nc-1;
                    active_imgviews(pp,kk:end) = 0;
                    if nc<2,
                        fprintf(1,'\nWarning: Only one active view left, Not able to generate 3D structure any more!\n');
                        break;
                    end;
                else
                    kth = (kk-1)*n_cam+pp;
                    area_mat(kth) = area_kk;
                    bounding_mat(:,kth) = boundbox_kk;
                    centroid_last(pp,:) = centroid_kk;
                    center2d_mat(:,kth) = centroid_kk';
                    foreground_cell{kth} = frame_kk(temp(2) : temp(4), temp(1) : temp(3));
                    % plot camera number
                    text(centroid_kk(1),centroid_kk(2),['\it\fontname{Arial}\fontsize{9}\color{green}' ...
                        num2str(pp)],'HorizontalAlignment','center');
                    rectangle('Position',boundbox_kk,'edgecolor','b','linewidth',1);
                end;
            end;
        end;
        if nc<2,
            break;
        end;
        title(['Check the foreground of frame ' num2str(kk) ':']);
        set(gca,'position',[0 0 1 1]);
        resolution = round(ny/height);   % dpi
        set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
        drawnow;
        if saveIM,
            save_name = [imgdir '/loc_' imgbase strnum_frame{kk} '.jpg'];
            print(2,'-djpeg',['-r' num2str(resolution)],save_name);
        end;
    end;
    fprintf(1,'\nDone with ROI extraction of all frames!\n');
    
else
    
    fprintf(1,'\nCompute image geometric feature under multiple-different-cameras configuration!\n');
    for pp = 1:n_cam,
        imgdir = imgdir_cell{pp};
        imgbase = imgbase_cell{pp};
        strnum_frame = strnum_frameNcam{pp};
        imPreNumf = [imgdir '/' fgprefix imgbase];
        active_view = active_imgviews(pp,:);
        ind_active = find(active_view);
        if isempty(ind_active),
            fprintf(1,'\nWarning: No frame is active in camera %d.\n',pp);
            continue;
        end;
        kk = ind_active(1);
        % load the 1st image of camera pp
        frame_kk = imread([imPreNumf strnum_frame{kk} '.png']);
        % geometric feature
        [y,x] = find(frame_kk);
        xx = [x,y]';
        area_kk = size(xx,2);
        temp = [min(xx,[],2); max(xx,[],2)];
        boundbox_kk = [temp(1:2)-0.5; temp(3:4)-temp(1:2)+1];
        distance_toler = max(boundbox_kk(3:4))*tol_factor;      % the longer edge of boundingbox
        centroid_kk = sum(xx,2)/area_kk;
        centroid_last = centroid_kk;
        threshold_area = area_kk*ndarea;
        kth = (kk-1)*n_cam+pp;
        bounding_mat(:,kth) = boundbox_kk;
        center2d_mat(:,kth) = centroid_kk;
        area_mat(kth) = area_kk;
        foreground_cell{kth} = frame_kk(temp(2) : temp(4), temp(1) : temp(3));
        
        figure(2); hold off;
        image(frame_kk);
        colormap(gray(2));      % colormap for binary image
        axis image;
        hold on;
        title(['Check the foreground of (camera ' num2str(pp) ', frame ' num2str(kk) '):']);
        plot(centroid_kk(1),centroid_kk(2),'g+');
        % plot boundingbox
        rectangle('Position',boundbox_kk,'edgecolor','b','linewidth',1);
        set(gca,'position',[0 0 1 1]);
        nx = imsize(1,pp);
        ny = imsize(2,pp);
        resolution = round(ny/height);
        set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
        drawnow;
        if saveIM,
            save_name = [imgdir '/loc_' imgbase strnum_frame{kk} '.jpg'];
            print(2,'-djpeg',['-r' num2str(resolution)],save_name);
        end;
        % load geometric feature of all image frames
        for kk = ind_active(2:end),
            frame_kk = imread([imPreNumf strnum_frame{kk} '.png']);
            % geometric feature
            [y,x] = find(frame_kk);
            xx = [x,y]';
            area_kk = size(xx,2);
            temp = [min(xx,[],2); max(xx,[],2)];
            boundbox_kk = [temp(1:2)-0.5; temp(3:4)-temp(1:2)+1];
            centroid_kk = sum(xx,2)/area_kk;
            figure(2); hold off;
            image(frame_kk);
            colormap(gray(2));      % colormap for binary image
            axis image;
            hold on;
            title(['check the foreground of (camera ' num2str(pp) ', frame ' num2str(kk) '):']);
            plot(centroid_kk(1),centroid_kk(2),'g+');
            % find the closest center of blobs to the clicked postion
            if  norm(centroid_kk-centroid_last)>distance_toler || area_kk < threshold_area,
                fprintf(1,'\nWarning: The foreground of (camera %d, frame %d) is abnormal!\n',pp,kk);
                active_view(kk) = 0;
            else
                kth = (kk-1)*n_cam+pp;
                bounding_mat(:,kth) = boundbox_kk;
                centroid_last = centroid_kk;
                center2d_mat(:,kth) = centroid_kk;
                area_mat(kth) = area_kk;
                foreground_cell{kth} = frame_kk(temp(2) : temp(4), temp(1) : temp(3));
                % plot boundingbox
                rectangle('Position',boundbox_kk,'edgecolor','b','linewidth',1);
            end;
            figure(2);
            set(gca,'position',[0 0 1 1]);
            set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
            drawnow;
            if saveIM,
                save_name = [imgdir '/loc_' imgbase strnum_frame{kk} '.jpg'];
                print(2,'-djpeg',['-r' num2str(resolution)],save_name);
            end;
        end;
        active_imgviews(pp,:) = active_view;
        fprintf(1,'\nDone with camera %d!\n',pp);
    end;
end;

active_images = sum(active_imgviews,1)>=2;
ind_active = find(active_images);
if isempty(ind_active),
    fprintf(1,'\nERROR: It takes at least two camera views to caculate 3d! No frame is available!\n');
    return;
end;

if ~exist('n_sublevel','var'),
    fprintf(1,'\nPlease input number of sublevels for visual hull refinement of level 1...\n');
    n_sublevel = input('Number of sublevels for level1: ([]=3) ');
    if isempty(n_sublevel),
        n_sublevel = 3;
    else
        n_sublevel = round(n_sublevel);
    end;
end;

if ~exist('n_subdiv','var') || ~isscalar(n_subdiv) || n_subdiv<2,
    fprintf(1,'\nThe subdivision number for sublevels should be set to 2 as recommended!\n');
    flag = 1;
    while flag,
        n_subdiv = input('Subdivision number for sublevels: ([]=2) ');
        if isempty(n_subdiv),
            n_subdiv = 2;
            flag = 0;
        else
            n_subdiv = round(n_subdiv);
            flag = (~isscalar(n_subdiv) || n_subdiv<2);
            if flag,
                fprintf(1,'\nThe subdivision number is less than 2! Please input again!\n');
            end;
        end;
    end;
end;

if ~exist('nfpu','var'),
    fprintf(1,'\nSome (mostly final) level will be divided into small units in case of memory deficiency!\n');
    nfpu = input('Set number of frames to save separately in one unit: ([]=50) ');
    if isempty(nfpu),
        nfpu = 50;
    end;
    n_unit = ceil(n_frame/nfpu);
end;

if ~exist('n_subdiv2','var') || ~isscalar(n_subdiv2) || n_subdiv2<2,
    fprintf(1,'\nSet subdivision number for the final level to compute 3D points cloud!\n');
    flag = 1;
    while flag,
        n_subdiv2 = input('Subdivision number to compute point cloud: ([]=10) ');
        if isempty(n_subdiv2),
            n_subdiv2 = 10;
            flag = 0;
        else
            n_subdiv2 = round(n_subdiv2);
            flag = (~isscalar(n_subdiv2) || n_subdiv2<2);
            if flag,
                fprintf(1,'\nThe subdivision number is less than 2! Please input again!\n');
            end;
        end;
    end;
end;

% compute level 0
compute_level0;

% compute level 1
compute_level1;

% saving environment variables
load_imgdir;
fprintf(1,'\nSaving environment variables as ''visualhull_environment.mat''.\n');
save_name = [imgdir '/visualhull_environment.mat'];
string_save = ['save ' save_name ' calib_mode fgprefix n_frame n_cam Tcc Omcc handcc ind_active active_images' ...
    ' active_imgviews character_len nxyz boundsize img_position maxbound bounding_mat area_mat foreground_cell' ...
    ' center2d_mat center3d_mat ex err_std ndigit cface height subnxyz level bricksize n_sublevel n_subdiv n_subdiv2 nfpu n_unit'];
if exist('Omcc_error','var'),
    string_save = [string_save ' Omcc_error Tcc_error'];
end;
if calib_mode,
    string_save = [string_save ' imgbase imgfmt strnum_frame nx ny fc cc alpha_c kc'];
    if exist('fc_error','var'),
        string_save = [string_save ' fc_error cc_error alpha_c_error kc_error'];
    end;
else
    string_save = [string_save ' imgbase_cell imgfmt_cell strnum_frameNcam' ...
        ' imsize fc_mat cc_mat alpha_vec kc_mat'];
    if exist('fc_mat_error','var'),
        string_save = [string_save ' fc_mat_error cc_mat_error alpha_vec_error kc_mat_error'];
    end;
end;
eval(string_save);

fprintf(1,'\nSaving the 1st level of visual hull data as ''visual_hull_level1.mat''.\n');
save_name = [imgdir '/visual_hull_level' num2str(level) '.mat'];
string_save = ['save ' save_name ' IN_hull ON_hull'];
eval(string_save);
fprintf(1,'done...\nNow refine the visual hull ...\n');

% finish the carving process
refine_visualhull;
