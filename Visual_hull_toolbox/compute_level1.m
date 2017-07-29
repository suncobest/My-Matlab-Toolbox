%% level 1
level =1;

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

if ~exist('subnxyz','var') || length(subnxyz)~=3,
    fprintf(1,['\nPlease input the subdivision number for the bounding box in the 1st level,\n' ...
        'it is recomended to have cuboid with equal side after subdivision ...\n']);
    flag = 1;
    while flag,
        subnxyz = input('subnxyz = ([]=[3;3;3])');
        if isempty(subnxyz),
            subnxyz = [3;3;3];
            flag=0;
        elseif length(subnxyz)==1,
            subnxyz = subnxyz(ones(3,1));
            flag=0;
        else
            subnxyz = subnxyz(:);
            flag = length(subnxyz)~=3;
            if flag,
                fprintf(1,'Unexpected input!');
            end;
        end;
    end;
end;
% Brick size of 1st level
bricksize = boundsize./subnxyz;
% number of bricks in the 1st level
nbricks = prod(subnxyz);
nverts = nbricks*8;

% voxel center of level 1 and their cube vertices
center = subvoxcenter(center3d_mat, subnxyz, boundsize);
vertices = gen_cuboids(center, bricksize);

IN_hull = cell(1, n_frame);
ON_hull = IN_hull;
if calib_mode,
    xx = NaN(2*n_cam, n_frame*nverts);
    for pp = 1:n_cam,
        om = Omcc(:,pp);
        T = Tcc(:,pp);
        hand = handcc(pp);
        x_kk = project_points_mirror2(vertices,om,T,hand,fc,cc,zeros(5,1),alpha_c);
        ii = (pp-1)*2;
        xx(ii+1:ii+2, :) = x_kk;
    end;
    
    if ~exist('imPreNumf','var'),
        imPreNumf = [imgdir '/' fgprefix imgbase];
    end;
    if ~exist('resolution','var'),
        resolution = round(ny/height);   % dpi
    end;
    if ~exist('tpx','var') || ~exist('tpy','var'),
        tpy = ny/20;
        tpx = nx-tpy; % nx/2;
    end;
    
    % plot projection of level 1 in every frame
    for kk = ind_active,
        jj = (kk-1)*nverts;
        % if a view is inactive, then all points cannot be eliminated in that view
        % default no cube in hull, all cubes on the border
        inner = true(n_cam, nbricks);
        border = inner;
        for pp = 1:n_cam,
            if active_imgviews(pp,kk),
                kth = (kk-1)*n_cam+pp;
                ii = (pp-1)*2;
                x_kk = xx(ii+1:ii+2, jj+1:jj+nverts)-repmat(bounding_mat(1:2,kth),1,nverts)+1.5;
                [in, on] = convhull_in_region(foreground_cell{kth}, x_kk, 8);
                inner(pp,:) = in;
                border(pp,:) = on | in;
            end;
        end;
        % visual hull data
        inner = all(inner,1);
        border = all(border,1) & (~inner);
        mm = (kk-1)*nbricks;
        X_kk = center(:,mm+1:mm+nbricks);
        IN_hull{kk} = X_kk(:, inner);
        ON_hull{kk} = X_kk(:, border);
        
        if saveIM,
            % face index of on border cuboids
            nb = sum(border);
            subface = reshape(repmat(cface(:),1,nb)+repmat(8*(0:nb-1),24,1), 4, []);
            % draw 2D projection of border cuboids in all views
            frame_kk = false(ny,nx);
            for pp = 1:n_cam,
                if active_imgviews(pp,kk),
                    kth = (kk-1)*n_cam+pp;
                    temp = bounding_mat(:,kth);
                    temp(1:2) = ceil(temp([2 1]));
                    frame_kk(temp(1) : temp(1)+temp(4)-1, temp(2) : temp(2)+temp(3)-1) = foreground_cell{kth};
                end;
            end;
            figure(2); hold off;
            image(frame_kk);
            colormap(gray(2));      % colormap for binary image
            hold on;
            framenb = sprintf(['%0' ndigit 'd'],kk);
            text(tpx,tpy,['\it\fontname{Arial}\fontsize{9}\color{white}Level \color{yellow}' ...
                num2str(level) '\color{white} - \color{yellow}',framenb],'HorizontalAlignment','right');
            xxj = xx(:, jj+1:jj+nverts)+1;
            xxj = xxj(:, reshape(border(ones(8,1),:),1,[]));          % logical indices
            for pp = 1:n_cam,
                if active_imgviews(pp,kk),
                    ii = (pp-1)*2;
                    x_kk = xxj(ii+1:ii+2, :);
                    patch('Faces', subface', 'Vertices', x_kk', 'FaceColor', 'b', 'EdgeColor','r','FaceAlpha',0.1);
                end;
            end;
            % set axis position in window, position = [left, bottom, width, height]
            set(gca,'position',[0 0 1 1]);
            set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
            drawnow;
            save_name = [imgdir '/level' num2str(level) '_' framenb '.jpg'];
            print(2,'-djpeg',['-r' num2str(resolution)],save_name);
            
            % plot visual hull and perspective image in 3D space
            figure(3); hold off;
            X_kk = [IN_hull{kk},ON_hull{kk}];
            plot3(X_kk(1,:),X_kk(3,:),-X_kk(2,:),'g.');
            hc = center3d_mat(:,kk);
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
            X_kk = X_kk(:, reshape(border(ones(8,1),:),1,[]));    % border boxes
            patch('Faces', subface', 'Vertices', [X_kk(1,:); X_kk(3,:); -X_kk(2,:)]', 'FaceColor', 'b','EdgeColor','r','FaceAlpha',0.2);
            axis equal tight vis3d off;
            view(az,el);
            %  lightangle(az-90,el);  set(3, 'renderer','zbuffer');
            set(3,'color',[0 0 0], 'PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
            drawnow;
            save_name = [imgdir '/3D_level' num2str(level) '_' framenb '.jpg'];
            print(3,'-djpeg',['-r' num2str(resolution)],save_name);
        end;
    end;
    
else
    
    xx = NaN(2*n_cam, n_frame*nverts);
    for pp = 1:n_cam,
        fc = fc_mat(:,pp);
        cc = cc_mat(:,pp);
        alpha_c = alpha_vec(pp);
        om = Omcc(:,pp);
        T = Tcc(:,pp);
        hand = handcc(pp);
        x_kk = project_points_mirror2(vertices,om,T,hand,fc,cc,zeros(5,1),alpha_c);
        ii = (pp-1)*2;
        xx(ii+1:ii+2, :) = x_kk;
    end;
    
    % plot projection of level 1 in every frame
    for kk = ind_active,
        jj = (kk-1)*nverts;
        inner = true(n_cam, nbricks);
        border = inner;
        for pp = 1:n_cam,
            if active_imgviews(pp,kk),
                kth = (kk-1)*n_cam+pp;
                ii = (pp-1)*2;
                x_kk = xx(ii+1:ii+2, jj+1:jj+nverts)-repmat(bounding_mat(1:2,kth),1,nverts)+1.5;
                [in, on] = convhull_in_region(foreground_cell{kth}, x_kk, 8);
                inner(pp,:) = in;
                border(pp,:) = on | in;
            end;
        end;
        inner = all(inner,1);
        border = all(border,1) & (~inner);
        mm = (kk-1)*nbricks;
        X_kk = center(:,mm+1:mm+nbricks);
        IN_hull{kk} = X_kk(:, inner);
        ON_hull{kk} = X_kk(:, border);
        
        if saveIM,
            framenb = sprintf(['%0' ndigit 'd'],kk);
            nb = sum(border);
            subface = reshape(repmat(cface(:),1,nb)+repmat(8*(0:nb-1),24,1), 4, []);
            xxj = xx(:, jj+1:jj+nverts)+1;
            xxj = xxj(:, reshape(border(ones(8,1),:),1,[]));          % logical indices
            % plot visual hull and perspective image in 3D space
            figure(3); hold off;
            X_kk = [IN_hull{kk},ON_hull{kk}];
            plot3(X_kk(1,:),X_kk(3,:),-X_kk(2,:),'g.');
            hc = center3d_mat(:,kk);
            for pp = 1:n_cam,
                if active_imgviews(pp,kk),
                    kth = (kk-1)*n_cam+pp;
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
                    text(tpx,tpy,['\it\fontname{Arial}\fontsize{9}\color{white}Level \color{yellow}' ...
                        num2str(level) '\color{white} - Cam \color{yellow}' num2str(pp) ...
                        '\color{white} - \color{yellow}', framenb],'HorizontalAlignment','right');
                    ii = (pp-1)*2;
                    x_kk = xxj(ii+1:ii+2, :);
                    patch('Faces', subface', 'Vertices', x_kk', 'FaceColor', 'b', 'EdgeColor','r','FaceAlpha',0.1);
                    % set axis position in window, position = [left, bottom, width, height]
                    set(gca,'position',[0 0 1 1]);
                    set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                    drawnow;
                    save_name = [imgdir_cell{pp} '/level' num2str(level) '_cam' num2str(pp) '_' framenb '.jpg'];
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
            X_kk = X_kk(:, reshape(border(ones(8,1),:),1,[]));
            patch('Faces', subface', 'Vertices', [X_kk(1,:); X_kk(3,:); -X_kk(2,:)]', 'FaceColor', 'b','EdgeColor','r','FaceAlpha',0.2);
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
            save_name = [imgdir '/3D_level' num2str(level) '_' framenb '.jpg'];
            print(3,'-djpeg',['-r' num2str(resolution)],save_name);
        end;
    end;
end;
