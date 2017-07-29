% refine visual hull cuboids in separate file units
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

% bricksize of last level
bricksize = (boundsize./subnxyz)/(n_subdiv^(level-1));

% initial 2*2 matrix for perspective image corner (surface function)
xc = [0 1;0 1];
yc = [0 0;1 1];
zc = [1 1;1 1];
% view axis for 3D space
az = 50;
el = 45;

for nl = level+1 : 1+n_sublevel,
    % Brick size of current level
    subbricksize = bricksize/n_subdiv;
    for count = 1:n_unit,
        base = (count-1)*nfpu;
        ind = base+1 : min(base+nfpu, n_frame);
        ind_active = find(active_images(ind));
        nc = sprintf(['%0' ndigit 'd'],count);
        % load part data of last level
        save_name = [imgdir '/level' num2str(nl-1) '_part' nc '.mat'];
        load(save_name);
        npts_inhull = cellfun(@(x) size(x,2), IN_hull);
        npts_onhull = cellfun(@(x) size(x,2), ON_hull);
        % number of bricks and verts of current level
        nbricks_inhull = npts_inhull*n_subdiv^3;
        nbricks_onhull = npts_onhull*n_subdiv^3;
        %     nverts_inhull = nbricks_inhull*8;
        nverts_onhull = nbricks_onhull*8;
        sum_nverts = sum(nverts_onhull);
        
        % voxel center of current level and their cube vertices
        INmat = cell2mat(IN_hull);
        ONmat = cell2mat(ON_hull);
        if isempty(ONmat),
            fprintf(1,'\nERROR: No border cuboids found in data %s!\n', save_name);
            return;
        end;
        center_inhull = [];
        if ~isempty(INmat),
            center_inhull = subvoxcenter(INmat, n_subdiv, bricksize);
            %         vertices_inhull = gen_cuboids(center_inhull, subbricksize);
        end;
        center_onhull = subvoxcenter(ONmat, n_subdiv, bricksize);
        vertices_onhull = gen_cuboids(center_onhull, subbricksize);
        % process inhull points
        IN_hull = cell(size(IN_hull));
        ON_hull = IN_hull;
        if ~isempty(center_inhull),
            jj2 = 0;
            for kk = ind_active,
                jj1 = jj2+1;
                jj2 = jj2 + nbricks_inhull(kk);
                IN_hull{kk} = center_inhull(:, jj1:jj2);
            end;
        end;
        
        if calib_mode,
            xx = NaN(2*n_cam, sum_nverts);
            for pp = 1:n_cam,
                om = Omcc(:,pp);
                T = Tcc(:,pp);
                hand = handcc(pp);
                x_kk = project_points_mirror2(vertices_onhull,om,T,hand,fc,cc,zeros(5,1),alpha_c);
                ii = (pp-1)*2;
                xx(ii+1:ii+2, :) = x_kk;
            end;
            
            % plot projection of 3D cuboids in every frame
            jj2 = 0; mm2 = 0;
            for kk = ind_active,
                bk = base+kk;
                jj1 = jj2+1;
                jj2 = jj2 + nverts_onhull(kk);
                % if a view is inactive, then all points cannot be eliminated in that view
                % default no cube in hull, all cubes on the border
                inner = true(n_cam, nbricks_onhull(kk));
                border = inner;
                for pp = 1:n_cam,
                    if active_imgviews(pp,bk),
                        kth = (bk-1)*n_cam+pp;
                        ii = (pp-1)*2;
                        x_kk = xx(ii+1:ii+2, jj1:jj2)-repmat(bounding_mat(1:2,kth),1,nverts_onhull(kk))+1.5;
                        [in, on] = convhull_in_region(foreground_cell{kth}, x_kk, 8);
                        inner(pp,:) = in;
                        border(pp,:) = on | in;
                    end;
                end;
                % visual hull data
                inner = all(inner,1);
                border = all(border,1) & (~inner);
                mm1 = mm2+1;
                mm2 = mm2 + nbricks_onhull(kk);
                X_kk = center_onhull(:, mm1:mm2);
                IN_hull{kk} = [IN_hull{kk}, X_kk(:, inner)];
                ON_hull{kk} = X_kk(:, border);
                
                if saveIM,
                    % face index of on border cuboids
                    nb = sum(border);
                    subface = reshape(repmat(cface(:),1,nb)+repmat(8*(0:nb-1),24,1), 4, []);
                    % draw 2D projection of border cuboids in all views
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
                    text(tpx,tpy,['\it\fontname{Arial}\fontsize{9}\color{white}Level \color{yellow}' ...
                        num2str(nl) '\color{white} - Frame \color{yellow}',framenb],'HorizontalAlignment','right');
                    xxj = xx(:, jj1:jj2)+1;
                    xxj = xxj(:, reshape(border(ones(8,1),:),1,[]));          % logical indices
                    for pp = 1:n_cam,
                        if active_imgviews(pp,bk),
                            ii = (pp-1)*2;
                            x_kk = xxj(ii+1:ii+2, :);
                            patch('Faces', subface', 'Vertices', x_kk', 'FaceColor', 'b', 'EdgeColor','r','FaceAlpha',0.05);
                        end;
                    end;
                    % set axis position in window, position = [left, bottom, width, height]
                    set(gca,'position',[0 0 1 1]);
                    set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                    drawnow;
                    save_name = [imgdir '/level' num2str(nl) '_frame' framenb '.jpg'];
                    print(2,'-djpeg',['-r' num2str(resolution)],save_name);
                    
                    % plot visual hull and perspective image in 3D space
                    figure(3); hold off;
                    X_kk = [IN_hull{kk},ON_hull{kk}];
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
                    X_kk = vertices_onhull(:, jj1:jj2);
                    X_kk = X_kk(:, reshape(border(ones(8,1),:),1,[]));    % border boxes
                    patch('Faces', subface', 'Vertices', [X_kk(1,:); X_kk(3,:); -X_kk(2,:)]', 'FaceColor', 'b','EdgeColor','r','FaceAlpha',0.1);
                    axis equal tight vis3d off;
                    view(az,el);
                    %  lightangle(az-90,el);  set(3, 'renderer','zbuffer');
                    set(3,'color',[0 0 0], 'PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                    drawnow;
                    save_name = [imgdir '/3D_level' num2str(nl) '_frame' framenb '.jpg'];
                    print(3,'-djpeg',['-r' num2str(resolution)],save_name);
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
                x_kk = project_points_mirror2(vertices_onhull,om,T,hand,fc,cc,zeros(5,1),alpha_c);
                ii = (pp-1)*2;
                xx(ii+1:ii+2, :) = x_kk;
            end;
            
            % plot projection of 3D cuboids in every frame
            jj2 = 0; mm2=0;
            for kk = ind_active,
                bk = base+kk;
                jj1 = jj2+1;
                jj2 = jj2 + nverts_onhull(kk);
                inner = true(n_cam, nbricks_onhull(kk));
                border = inner;
                for pp = 1:n_cam,
                    if active_imgviews(pp,bk),
                        kth = (bk-1)*n_cam+pp;
                        ii = (pp-1)*2;
                        x_kk = xx(ii+1:ii+2, jj1:jj2)-repmat(bounding_mat(1:2,kth),1,nverts_onhull(kk))+1.5;
                        [in, on] = convhull_in_region(foreground_cell{kth}, x_kk, 8);
                        inner(pp,:) = in;
                        border(pp,:) = on | in;
                    end;
                end;
                inner = all(inner,1);
                border = all(border,1) & (~inner);
                mm1 = mm2+1;
                mm2 = mm2 + nbricks_onhull(kk);
                X_kk = center_onhull(:, mm1:mm2);
                IN_hull{kk} = [IN_hull{kk},X_kk(:, inner)];
                ON_hull{kk} = X_kk(:, border);
                
                if saveIM,
                    framenb = sprintf(['%0' ndigit 'd'],bk);
                    nb = sum(border);
                    subface = reshape(repmat(cface(:),1,nb)+repmat(8*(0:nb-1),24,1), 4, []);
                    xxj = xx(:, jj1:jj2)+1;
                    xxj = xxj(:, reshape(border(ones(8,1),:),1,[]));          % logical indices
                    % plot visual hull and perspective image in 3D space
                    figure(3); hold off;
                    X_kk = [IN_hull{kk},ON_hull{kk}];
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
                            text(tpx,tpy,['\it\fontname{Arial}\fontsize{9}\color{white}Level \color{yellow}' ...
                                num2str(nl) '\color{white} - Camera \color{yellow}' num2str(pp) ...
                                '\color{white} - Frame \color{yellow}', framenb],'HorizontalAlignment','right');
                            ii = (pp-1)*2;
                            x_kk = xxj(ii+1:ii+2, :);
                            patch('Faces', subface', 'Vertices', x_kk', 'FaceColor', 'b', 'EdgeColor','r','FaceAlpha',0.05);
                            % set axis position in window, position = [left, bottom, width, height]
                            set(gca,'position',[0 0 1 1]);
                            set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                            drawnow;
                            save_name = [imgdir_cell{pp} '/level' num2str(nl) '_camera' num2str(pp) '_frame' framenb '.jpg'];
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
                    X_kk = vertices_onhull(:, jj1:jj2);
                    X_kk = X_kk(:, reshape(border(ones(8,1),:),1,[]));
                    patch('Faces', subface', 'Vertices', [X_kk(1,:); X_kk(3,:); -X_kk(2,:)]', 'FaceColor', 'b','EdgeColor','r','FaceAlpha',0.1);
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
                    save_name = [imgdir '/3D_level' num2str(nl) '_frame' framenb '.jpg'];
                    print(3,'-djpeg',['-r' num2str(resolution)],save_name);
                end;
            end;
            imgdir = imgdir_cell{1};
        end;
        fprintf(1,'\nSaving parted visualhull data ''level%d_part%s.mat''.\n', nl,nc);
        save_name = [imgdir '/level' num2str(nl) '_part' nc '.mat'];
        eval(['save ' save_name ' IN_hull ON_hull']);
    end;
    bricksize = subbricksize;
    level = level+1;
end;
ind_active = find(active_images);
