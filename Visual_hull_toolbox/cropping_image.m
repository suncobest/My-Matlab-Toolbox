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
    fprintf(1,'\nNo file ''load_imgdir.m'' found in the code directory!\nPlease compute visual hull first!\n');
    return;
end;

save_name = [imgdir '/kinematic_animation.mat'];
if exist(save_name,'file')==2,
    fprintf(1,'\nLoading kinematics parameters file ''kinematic_animation.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\nERROR: ''kinematic_animation.mat'' not found!\n');
    return;
end;

save_name = [imgdir '/visualhull_environment.mat'];
if exist(save_name, 'file')==2,
    fprintf(1,'\nLoading environment variable file ''visualhull_environment.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\nERROR: ''visualhull_environment.mat'' not found!\n');
    return;
end;

a = min(member_center(:,1,ind_active),[],3)-maxbound/2;
b = max(member_center(:,1,ind_active),[],3)+maxbound/2;
center = (a+b)/2;     % center of 3D space
ax = max(b-a)/2;
img_position =  max(ax, img_position)-maxbound/10;
for pp = 1:n_cam,
    om = Omcc(:,pp);
    T = Tcc(:,pp);
    hand = handcc(pp);
    if ~calib_mode,    % multiple cameras
        fc = fc_mat(:,pp);
        cc = cc_mat(:,pp);
        alpha_c = alpha_vec(pp);
    end;
    % direction of camera frame axes in world frame: Rkk
    % camera center in world frame: -Rkk*T
    Rkk = rodrigues(-om);       % transpose
    if hand~=1,
        Rkk(3,:) = -Rkk(3,:);
    end;
    % positon of pespective image center
    vect = center+Rkk*T;
    XX = center+img_position*vect/norm(vect);
    % position of pespective image corners (o, y, xy, x)
    XX = XX(:,ones(1,4))+(Rkk(:,1)*[-1,-1,1,1]+Rkk(:,2)*[-1,1,1,-1])*ax;
    % project the four corner on image
    xx = project_points_mirror2(XX,om,T,hand,fc,cc,zeros(5,1),alpha_c)+1;
    % transform to matlab 3D axes for surface function (o, y, x, xy)
    if calib_mode,
        for kk=ind_active,
            if active_imgviews(pp,kk),
                framenb = strnum_frame{kk};
                % load original images
                frame_kk = imread([imgdir '/' imgbase framenb '.' imgfmt]);
                frame_kk = homography_image(frame_kk, xx);
                imwrite(frame_kk,[imgdir '/cam' num2str(pp) '_old_'   framenb '.jpg'],'jpg');
                % load images with reprojection
                frame_kk = imread([imgdir '/projection_' imgbase framenb '.jpg']);
                frame_kk = homography_image(frame_kk, xx);
                imwrite(frame_kk,[imgdir '/cam' num2str(pp) '_prj_'  framenb '.jpg'],'jpg');
                % load B&W foregound images
                frame_kk = false(ny,nx);
                kth = (kk-1)*n_cam+pp;
                temp = bounding_mat(:,kth);
                temp(1:2) = ceil(temp([2 1]));
                frame_kk(temp(1) : temp(1)+temp(4)-1, temp(2) : temp(2)+temp(3)-1) = foreground_cell{kth};
                frame_kk = homography_image(frame_kk, xx);
                imwrite(frame_kk,[imgdir '/cam' num2str(pp) '_bw_' framenb '.png'],'png');
            end;
        end;
    else
        imgdir = imgdir_cell{pp};
        imgbase = imgbase_cell{pp};
        imgfmt = imgfmt_cell{pp};
        nx = imsize(1,pp);
        ny = imsize(2,pp);
        for kk=ind_active,
            if active_imgviews(pp,kk),
                framenb = strnum_frameNcam{pp}{kk};
                % load original images
                frame_kk = imread([imgdir '/' imgbase framenb '.' imgfmt]);
                frame_kk = homography_image(frame_kk, xx);
                imwrite(frame_kk,[imgdir '/cam' num2str(pp) '_old_'   framenb '.jpg'],'jpg');
                % load images with reprojection
                frame_kk = imread([imgdir '/projection_' imgbase framenb '.jpg']);
                frame_kk = homography_image(frame_kk, xx);
                imwrite(frame_kk,[imgdir '/cam' num2str(pp) '_prj_'  framenb '.jpg'],'jpg');
                % load B&W foregound images
                frame_kk = false(ny,nx);
                kth = (kk-1)*n_cam+pp;
                temp = bounding_mat(:,kth);
                temp(1:2) = ceil(temp([2 1]));
                frame_kk(temp(1) : temp(1)+temp(4)-1, temp(2) : temp(2)+temp(3)-1) = foreground_cell{kth};
                frame_kk = homography_image(frame_kk, xx);
                imwrite(frame_kk,[imgdir '/cam' num2str(pp) '_bw_' framenb '.png'],'png');
            end;
        end;
    end;
end;
