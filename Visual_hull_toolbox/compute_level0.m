%% level 0
level = 0;
% vertex index of every cube face
cface = [1, 5, 1, 4, 2, 1;
    4, 6, 2, 8, 3, 5;
    3, 7, 6, 7, 7, 8;
    2, 8, 5, 3, 6, 4];
% string form of image number digit
ndigit = num2str(floor(log10(n_frame))+1);

if ~exist('img_position','var'),
    fprintf(1,['\nWhere do you want to put up the perspective images in 3D space?\n' ...
        'You can set the image back behind the visual hull in unit of characteristic length ...\n']);
    img_position = input('img_position = character_len * ([]=5) ');
    if isempty(img_position),
        img_position = 5*character_len;
    else
        img_position = img_position*character_len;
    end;
end;

if ~exist('MaxIter','var'),
    MaxIter = 10;       % number of iteration to compute 3d pints
end;

if ~exist('n_view','var'),
    n_view = n_frame*n_cam;
end;

if ~exist('saveIM','var'),
    saveIM = input('Draw and save 2D, 3D images or not? ([]=no, other=yes) ','s');
    saveIM = ~isempty(saveIM);
end;

% initial 2*2 matrix for perspective image corner (surface function)
xc = [0 1;0 1];
yc = [0 0;1 1];
zc = [1 1;1 1];
% view axis for 3D space
az = 50;
el = 45;

% number of vertices
nverts = 8;
% length of boundsize edge
maxbound = max(boundsize);
% rearrange center2d_mat to meet the requirement of function compute_structure
% dimension of (2, n_frame, n_cam), origin point is (0,0) instead of (1,1)
center = permute(reshape(center2d_mat, 2, n_cam, n_frame), [1, 3, 2])-1;

if calib_mode,
    center3d_mat = compute_structure2(center,Omcc,Tcc,handcc,fc,cc,zeros(5,1),alpha_c,MaxIter);
    % center of all 3D bounding boxes and their vertices
    vertices = gen_cuboids(center3d_mat,boundsize);
    X_kk = [center3d_mat, vertices];
    ex = NaN(2*n_cam, n_frame);
    xx = NaN(2*n_cam, nverts*n_frame);
    for pp = 1:n_cam,
        om = Omcc(:,pp);
        T = Tcc(:,pp);
        hand = handcc(pp);
        x_kk = project_points_mirror2(X_kk,om,T,hand,fc,cc,zeros(5,1),alpha_c);
        ii = (pp-1)*2;
        ex(ii+1:ii+2, :) = center(:,:,pp)-x_kk(:,1:n_frame);
        xx(ii+1:ii+2, :) = x_kk(:,n_frame+1:end);        % projection of cuboid vertices
    end;
    ex = reshape(ex,2,n_view);
    ex = ex(:,~isnan(ex(1,:)));               % get rid of inactive views;
    err_std = std(ex,0,2);
    
    if saveIM,
        if ~exist('imPreNumf','var'),
            imPreNumf = [imgdir '/' fgprefix imgbase];
        end;
        if ~exist('resolution','var'),
            resolution = round(ny/height);   % dpi
        end;
        % plot visual hull of level 0 in every frame
        % text position. unit: points
        tpy = ny/20;
        tpx = nx-tpy; % nx/2;
        for kk = ind_active,
            frame_kk = imread([imPreNumf strnum_frame{kk} '.png']);
            figure(2); hold off;
            image(frame_kk);
            colormap(gray(2));      % colormap for binary image
            hold on;
            framenb = sprintf(['%0' ndigit 'd'],kk);
            text(tpx,tpy,['\it\fontname{Arial}\fontsize{9}\color{white}Level \color{yellow}' ...
                num2str(level) '\color{white} - \color{yellow}',framenb],'HorizontalAlignment','right');
            jj = (kk-1)*nverts;
            for pp = 1:n_cam,
                if active_imgviews(pp,kk),
                    ii = (pp-1)*2;
                    x_kk = xx(ii+1:ii+2, jj+1:jj+nverts);
                    figure(2);
                    patch('Faces', cface', 'Vertices', x_kk'+1, 'FaceColor', 'b', 'EdgeColor','r','FaceAlpha',0.1);
                end;
            end;
            % set axis position in window, position = [left, bottom, width, height]
            figure(2);
            set(gca,'position',[0 0 1 1]);
            set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
            drawnow;
            save_name = [imgdir '/level' num2str(level) '_' framenb '.jpg'];
            print(2,'-djpeg',['-r' num2str(resolution)],save_name);
            
            % plot visual hull and perspective image in 3D space
            figure(3); hold off;
            hc = center3d_mat(:,kk);
            plot3(hc(1),hc(3),-hc(2),'g.');
            frame_kk = imread(save_name);
            for pp = 1:n_cam,
                if active_imgviews(pp,kk),
                    om = Omcc(:,pp);
                    T = Tcc(:,pp);
                    hand = handcc(pp);
                    % direction of camera frame axes in world frame: Rkk
                    % camera center in world frame: -Rkk*T
                    Rkk = rodrigues(-om);       % transpose
                    if hand~=1,
                        Rkk(3,:) = -Rkk(3,:);
                    end;
                    % positon of pespective image center
                    vect = hc+Rkk*T;
                    X_kk = hc+img_position*vect/norm(vect);
                    % position of pespective image corners (o, y, xy, x)
                    X_kk = X_kk(:,ones(1,4))+(Rkk(:,1)*[-1,-1,1,1]+Rkk(:,2)*[-1,1,1,-1])*maxbound*0.6;
                    % project the four corner on image
                    x_kk = project_points_mirror2(X_kk,om,T,hand,fc,cc,zeros(5,1),alpha_c);
                    Ikk = homography_image(frame_kk, x_kk+1, [], 0);   % image background is black
                    % transform to matlab 3D axes for surface function (o, y, x, xy)
                    X_kk = X_kk(:, [1,2,4,3]);
                    xc(:) = X_kk(1,:);
                    yc(:) = X_kk(2,:);
                    zc(:) = X_kk(3,:);
                    surface(xc, zc, -yc, Ikk,'FaceColor','texturemap','EdgeColor','none');
                end;
            end;
            % set 3D plot
            figure(3);
            X_kk = vertices(:, jj+1:jj+nverts);
            patch('Faces', cface', 'Vertices', [X_kk(1,:); X_kk(3,:); -X_kk(2,:)]', 'FaceColor', 'b','EdgeColor','r','FaceAlpha',0.2);
            axis equal tight vis3d off;
            view(az,el);
            %   lightangle(az-90,el);  set(3, 'renderer','zbuffer');
            set(3,'color',[0 0 0], 'PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
            drawnow;
            save_name = [imgdir '/3D_level' num2str(level) '_' framenb '.jpg'];
            print(3,'-djpeg',['-r' num2str(resolution)],save_name);
        end;
    end;
    
else
    
    center3d_mat = compute_structure2(center,Omcc,Tcc,handcc,fc_mat,cc_mat,zeros(5,n_cam),alpha_vec,MaxIter);
    % center of all 3D bounding boxes and their vertices
    vertices = gen_cuboids(center3d_mat,boundsize);
    X_kk = [center3d_mat, vertices];
    ex = NaN(2*n_cam, n_frame);
    xx = NaN(2*n_cam, nverts*n_frame);
    for pp = 1:n_cam,
        fc = fc_mat(:,pp);
        cc = cc_mat(:,pp);
        alpha_c = alpha_vec(pp);
        om = Omcc(:,pp);
        T = Tcc(:,pp);
        hand = handcc(pp);
        x_kk = project_points_mirror2(X_kk,om,T,hand,fc,cc,zeros(5,1),alpha_c);
        ii = (pp-1)*2;
        ex(ii+1:ii+2, :) = center(:,:,pp)-x_kk(:,1:n_frame);
        xx(ii+1:ii+2, :) = x_kk(:,n_frame+1:end);        % projection of cuboid vertices
    end;
    ex = reshape(ex,2,n_view);
    ex = ex(:,~isnan(ex(1,:)));               % get rid of inactive views;
    err_std = std(ex,0,2);
    
    if saveIM,
        % plot projection of level 0 in every frame
        for kk = ind_active,
            framenb = sprintf(['%0' ndigit 'd'],kk);
            jj = (kk-1)*nverts;
            % plot visual hull and perspective image in 3D space
            figure(3); hold off;
            hc = center3d_mat(:,kk);
            plot3(hc(1),hc(3),-hc(2),'g.');
            for pp = 1:n_cam,
                if active_imgviews(pp,kk),
                    imgdir = imgdir_cell{pp};
                    imgbase = imgbase_cell{pp};
                    strnum_frame = strnum_frameNcam{pp};
                    frame_kk = imread([imgdir '/' fgprefix imgbase strnum_frame{kk} '.png']);
                    % draw 2D projection of boundingbox in all camera views
                    figure(2); hold off;
                    image(frame_kk);
                    colormap(gray(2));      % colormap for binary image
                    hold on;
                    nx = imsize(1,pp);
                    ny = imsize(2,pp);
                    % text position. unit: points
                    tpy = ny/20;
                    tpx = nx-tpy;   % nx/2;
                    resolution = round(ny/height);   % dpi
                    text(tpx,tpy,['\it\fontname{Arial}\fontsize{9}\color{white}Level \color{yellow}' ...
                        num2str(level) '\color{white} - Cam \color{yellow}' num2str(pp) ...
                        '\color{white} - \color{yellow}',framenb],'HorizontalAlignment','right');
                    ii = (pp-1)*2;
                    x_kk = xx(ii+1:ii+2, jj+1:jj+nverts);
                    patch('Faces', cface', 'Vertices', x_kk'+1, 'FaceColor', 'b', 'EdgeColor','r','FaceAlpha',0.1);
                    % set axis position in window, position = [left, bottom, width, height]
                    set(gca,'position',[0 0 1 1]);
                    set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                    drawnow;
                    save_name = [imgdir '/level' num2str(level) '_cam' num2str(pp) '_' framenb '.jpg'];
                    print(2,'-djpeg',['-r' num2str(resolution)],save_name);
                    
                    % draw 3D
                    frame_kk = imread(save_name);
                    fc = fc_mat(:,pp);
                    cc = cc_mat(:,pp);
                    alpha_c = alpha_vec(pp);
                    om = Omcc(:,pp);
                    T = Tcc(:,pp);
                    hand = handcc(pp);
                    % direction of camera frame axes in world frame: Rkk
                    % camera center in world frame: -Rkk*T
                    Rkk = rodrigues(-om);       % transpose
                    if hand~=1,
                        Rkk(3,:) = -Rkk(3,:);
                    end;
                    % positon of pespective image center
                    vect = hc+Rkk*T;
                    X_kk = hc+img_position*vect/norm(vect);
                    % position of pespective image corners (o, y, xy, x)
                    X_kk = X_kk(:,ones(1,4))+(Rkk(:,1)*[-1,-1,1,1]+Rkk(:,2)*[-1,1,1,-1])*maxbound*0.6;
                    % project the four corner on image
                    x_kk = project_points_mirror2(X_kk,om,T,hand,fc,cc,zeros(5,1),alpha_c);
                    Ikk = homography_image(frame_kk, x_kk+1, [], 0);   % image background is black
                    % transform to matlab 3D axes for surface function (o, y, x, xy)
                    X_kk = X_kk(:, [1,2,4,3]);
                    xc(:) = X_kk(1,:);
                    yc(:) = X_kk(2,:);
                    zc(:) = X_kk(3,:);
                    figure(3);
                    surface(xc, zc, -yc, Ikk,'FaceColor','texturemap','EdgeColor','none');
                end;
            end;
            % set 3D plot
            figure(3);
            X_kk = vertices(:, jj+1:jj+nverts);
            patch('Faces', cface', 'Vertices', [X_kk(1,:); X_kk(3,:); -X_kk(2,:)]', 'FaceColor', 'b','EdgeColor','r','FaceAlpha',0.2);
            axis equal tight vis3d off;
            view(az,el);
            nx = imsize(1,1);
            ny = imsize(2,1);
            resolution = round(ny/height);
            %  lightangle(az-90,el);  set(3, 'renderer','zbuffer');
            set(3,'color',[0 0 0], 'PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
            drawnow;
            imgdir = imgdir_cell{1};
            save_name = [imgdir '/3D_level' num2str(level) '_' framenb '.jpg'];
            print(3,'-djpeg',['-r' num2str(resolution)],save_name);
        end;
    end;
end;
