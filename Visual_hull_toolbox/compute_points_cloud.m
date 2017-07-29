%%   compute points cloud
if ~exist('saveIM','var'),
    saveIM = input('Draw and save 2D, 3D images or not? ([]=no, other=yes) ','s');
    saveIM = ~isempty(saveIM);
end;

if ~exist('level','var') || ~exist('n_subdiv','var'),
    load_imgdir;
    save_name = [imgdir '/visualhull_environment.mat'];
    if exist(save_name, 'file')==2,
        fprintf(1,'\nLoading environment variable file ''visualhull_environment.mat'' ...\n');
        load(save_name);
    else
        fprintf(1,'\nERROR: ''visualhull_environment.mat'' not found!\n');
        return;
    end;
end;

% bricksize of final level
bricksize = (boundsize./subnxyz)/(n_subdiv^n_sublevel);

% initial 2*2 matrix for perspective image corner (surface function)
xc = [0 1;0 1];
yc = [0 0;1 1];
zc = [1 1;1 1];
% view axis for 3D space
az = 50;
el = 45;

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

if calib_mode,
    if ~exist('resolution','var'),
        resolution = round(ny/height);   % dpi
    end;
    if ~exist('tpx','var') || ~exist('tpy','var'),
        tpy = ny/20;
        tpx = nx-tpy; % nx/2;
    end;
end;

for count = 1:n_unit,
    base = (count-1)*nfpu;
    ind = base+1 : min(base+nfpu, n_frame);
    ind_active = find(active_images(ind));
    nc = sprintf(['%0' ndigit 'd'],count);
    % load part data of last level
    save_name = [imgdir '/level' num2str(level) '_part' nc '.mat'];
    load(save_name);
    npts_inhull = cellfun(@(x) size(x,2), IN_hull);
    npts_onhull = cellfun(@(x) size(x,2), ON_hull);
    % number of bricks and verts of current level
    nbricks_inhull = npts_inhull*n_subdiv2^3;
    nbricks_onhull = npts_onhull*n_subdiv2^3;
    sum_nverts = sum(nbricks_onhull);
    
    % voxel center of current level and their cube vertices
    INmat = cell2mat(IN_hull);
    ONmat = cell2mat(ON_hull);
    if isempty(ONmat),
        fprintf(1,'\nERROR: No border cuboids found in data %s!\n', save_name);
        return;
    end;
    center_inhull = [];
    if ~isempty(INmat),
        center_inhull = subvoxcenter(INmat, n_subdiv2, bricksize);
    end;
    center_onhull = subvoxcenter(ONmat, n_subdiv2, bricksize);
    % process inhull points
    X_cell = cell(size(IN_hull));
    if ~isempty(center_inhull),
        jj2 = 0;
        for kk = ind_active,
            jj1 = jj2+1;
            jj2 = jj2 + nbricks_inhull(kk);
            X_cell{kk} = center_inhull(:, jj1:jj2);
        end;
    end;
    
    if calib_mode,
        xx = NaN(2*n_cam, sum_nverts);
        for pp = 1:n_cam,
            om = Omcc(:,pp);
            T = Tcc(:,pp);
            hand = handcc(pp);
            x_kk = project_points_mirror2(center_onhull,om,T,hand,fc,cc,zeros(5,1),alpha_c);
            ii = (pp-1)*2;
            xx(ii+1:ii+2, :) = x_kk;
        end;
        
        % plot projection of 3D points in every frame
        jj2 = 0;
        for kk = ind_active,
            bk = base+kk;
            jj1 = jj2+1;
            jj2 = jj2 + nbricks_onhull(kk);
            % xx stores points in cuboid on border
            inner = true(n_cam, nbricks_onhull(kk));
            for pp = 1:n_cam,
                if active_imgviews(pp,bk),
                    kth = (bk-1)*n_cam+pp;
                    ii = (pp-1)*2;
                    x_kk = xx(ii+1:ii+2, jj1:jj2)-repmat(bounding_mat(1:2,kth),1,nbricks_onhull(kk))+1.5;
                    [in, on] = convhull_in_region(foreground_cell{kth}, x_kk, 1);
                    inner(pp,:) = on | in;
                end;
            end;
            % visual hull data
            inner = all(inner,1);
            X_kk = center_onhull(:, jj1:jj2);
            X_kk = X_kk(:, inner);
            X_cell{kk} = [X_cell{kk}, X_kk];

            if saveIM,
                % draw 2D projection of border points in all views
                frame_kk = false(ny,nx);
                for pp = 1:n_cam,
                    if active_imgviews(pp,bk),
                        kth = (bk-1)*n_cam+pp;
                        temp = bounding_mat(:,kth);
                        temp(1:2) = ceil(temp([2 1]));
                        frame_kk(temp(1) : temp(1)+temp(4)-1, temp(2) : temp(2)+temp(3)-1) = foreground_cell{kth};
                    end;
                end;
                figure(2); hold off;
                image(frame_kk);
                colormap(gray(2));      % colormap for binary image
                hold on;
                framenb = sprintf(['%0' ndigit 'd'],bk);
                text(tpx,tpy,['\it\fontname{Arial}\fontsize{9}\color{white}Points - ' ...
                    '\color{yellow}',framenb],'HorizontalAlignment','right');
                xxj = xx(:, jj1:jj2)+1;
                xxj = xxj(:,inner);         % logical indices
                for pp = 1:n_cam,
                    if active_imgviews(pp,bk),
                        ii = (pp-1)*2;
                        x_kk = xxj(ii+1:ii+2, :);
                        plot(x_kk(1,:),x_kk(2,:),'g.');
                    end;
                end;
                % set axis position in window, position = [left, bottom, width, height]
                set(gca,'position',[0 0 1 1]);
                set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                drawnow;
                save_name = [imgdir '/points_' framenb '.jpg'];
                print(2,'-djpeg',['-r' num2str(resolution)],save_name);
                
                % plot 3D points on border and perspective image in 3D space
                figure(3); hold off;
                plot3(X_kk(1,:),X_kk(3,:),-X_kk(2,:),'g.');
                hc = center3d_mat(:,bk);
                frame_kk = imread(save_name);
                for pp = 1:n_cam,
                    if active_imgviews(pp,bk),
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
                axis equal tight vis3d off;
                view(az,el);
                %  lightangle(az-90,el);  set(3, 'renderer','zbuffer');
                set(3,'color',[0 0 0], 'PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                drawnow;
                save_name = [imgdir '/3D_points_' framenb '.jpg'];
                print(3,'-djpeg',['-r' num2str(resolution)],save_name);
                saveas(3, [imgdir '/3D_points_' framenb],'fig');
            end;
        end;
        
    else
        
        xx = NaN(2*n_cam, sum_nverts);
        for pp = 1:n_cam,
            fc = fc_mat(:,pp);
            cc = cc_mat(:,pp);
            alpha_c = alpha_vec(pp);
            om = Omcc(:,pp);
            T = Tcc(:,pp);
            hand = handcc(pp);
            x_kk = project_points_mirror2(center_onhull,om,T,hand,fc,cc,zeros(5,1),alpha_c);
            ii = (pp-1)*2;
            xx(ii+1:ii+2, :) = x_kk;
        end;
        
        % plot projection of 3D points in every frame
        jj2 = 0;
        for kk = ind_active,
            bk = base+kk;
            jj1 = jj2+1;
            jj2 = jj2 + nbricks_onhull(kk);
            inner = true(n_cam, nbricks_onhull(kk));
            for pp = 1:n_cam,
                if active_imgviews(pp,bk),
                    kth = (bk-1)*n_cam+pp;
                    ii = (pp-1)*2;
                    x_kk = xx(ii+1:ii+2, jj1:jj2)-repmat(bounding_mat(1:2,kth),1,nbricks_onhull(kk))+1.5;
                    [in, on] = convhull_in_region(foreground_cell{kth}, x_kk, 1);
                    inner(pp,:) = on | in;
                end;
            end;
            inner = all(inner,1);
            X_kk = center_onhull(:, jj1:jj2);
            X_kk = X_kk(:,inner);
            X_cell{kk} = [X_cell{kk}, X_kk];

            if saveIM,
                framenb = sprintf(['%0' ndigit 'd'],bk);
                xxj = xx(:, jj1:jj2)+1;
                xxj = xxj(:, inner);          % logical indices
                % plot visual hull and perspective image in 3D space
                figure(3); hold off;
                plot3(X_kk(1,:),X_kk(3,:),-X_kk(2,:),'g.');
                hc = center3d_mat(:,bk);
                for pp = 1:n_cam,
                    if active_imgviews(pp,bk),
                        kth = (bk-1)*n_cam+pp;
                        nx = imsize(1,pp);
                        ny = imsize(2,pp);
                        frame_kk = false(ny,nx);
                        temp = bounding_mat(:,kth);
                        temp(1:2) = ceil(temp([2 1]));
                        frame_kk(temp(1) : temp(1)+temp(4)-1, temp(2) : temp(2)+temp(3)-1) = foreground_cell{kth};
                        % draw 2D projection of border cuboids in all camera views
                        figure(2); hold off;
                        image(frame_kk);
                        colormap(gray(2));      % colormap for binary image
                        hold on;
                        % text position. unit: points
                        tpy = ny/20;
                        tpx = nx-tpy;  % nx/2;
                        resolution = round(ny/height);   % dpi
                        text(tpx,tpy,['\it\fontname{Arial}\fontsize{9}\color{white}Points - Cam \color{yellow}' ...
                            num2str(pp) '\color{white} - \color{yellow}', framenb],'HorizontalAlignment','right');
                        ii = (pp-1)*2;
                        x_kk = xxj(ii+1:ii+2, :);
                        plot(x_kk(1,:),x_kk(2,:),'g.');
                        % set axis position in window, position = [left, bottom, width, height]
                        set(gca,'position',[0 0 1 1]);
                        set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                        drawnow;
                        save_name = [imgdir_cell{pp} '/points_cam' num2str(pp) '_' framenb '.jpg'];
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
                axis equal tight vis3d off;
                view(az,el);
                nx = imsize(1,1);
                ny = imsize(2,1);
                resolution = round(ny/height);
                %  lightangle(az-90,el);  set(3, 'renderer','zbuffer');
                set(3,'color',[0 0 0], 'PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                drawnow;
                % set the saving directory to the first camera view folder
                imgdir = imgdir_cell{1};
                save_name = [imgdir '/3D_points_' framenb '.jpg'];
                print(3,'-djpeg',['-r' num2str(resolution)],save_name);
                saveas(3, [imgdir '/3D_points_' framenb],'fig');
            end;
        end;
         imgdir = imgdir_cell{1};
    end;
    fprintf(1,'\nSaving parted 3D points cloud data ''3D_points_part%s.mat''.\n', nc);
    save_name = [imgdir '/3D_points_part' nc '.mat'];
    eval(['save ' save_name ' X_cell']);
end;

% points interval
bricksize = bricksize/n_subdiv2;
ind_active = find(active_images);

% saving environment variables
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
fprintf(1,'\nEnvironment variables is updated in ''visualhull_environment.mat''.\n');

if saveIM,
    fprintf(1,'\nRun script ''open([imgdir ''/3D_points_#kk.fig'']);'' to open a saved figure ...\n');
end;
