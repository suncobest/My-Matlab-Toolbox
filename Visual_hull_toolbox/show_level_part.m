if ~exist('SW','var'),
    % show points on hull
    SW = 1;
end;

if ~exist('draw_faces','var'),
    % show cuboids' faces
    draw_faces= 1;
end;

if ~exist('Rd','var'),
    % render 3D faces
    Rd= 1;
end;

if ~exist('pim','var'),
    % place 2D images in 3D space
    pim= 1;
end;

for count = 1:n_unit,
    base = (count-1)*nfpu;
    ind = base+1 : min(base+nfpu, n_frame);
    ind_active = find(active_images(ind));
    nc = sprintf(['%0' ndigit 'd'],count);
    % load part data of last level
    save_name = [imgdir '/level' num2str(Nb) '_part' nc '.mat'];
    load(save_name);
    
    if SW,
        npts = cellfun(@(x) size(x,2), ON_hull);
        XX = cell2mat(ON_hull);
    else
        npts = cellfun(@(x) size(x,2), IN_hull);
        XX = cell2mat(IN_hull);
    end;
    if isempty(XX),
        fprintf(1,'\nNo cuboids found!\n');
        return;
    end;
    vertices = gen_cuboids(XX, bricksize);
    nverts = npts*8;
    sum_npts = size(XX,2);
    
    if calib_mode,
        xx = NaN(2*n_cam, sum_npts);
        verts2d = NaN(2*n_cam, sum_npts*8);
        for pp = 1:n_cam,
            om = Omcc(:,pp);
            T = Tcc(:,pp);
            hand = handcc(pp);
            x_kk = project_points_mirror2([XX,vertices],om,T,hand,fc,cc,zeros(5,1),alpha_c);
            ii = (pp-1)*2;
            xx(ii+1:ii+2, :) = x_kk(:,1:sum_npts);
            verts2d(ii+1:ii+2, :) = x_kk(:,sum_npts+1:end);
        end;
        
        jj2 = 0; mm2 = 0;
        for kk = ind_active,
            bk = base+kk;
            framenb = sprintf(['%0' ndigit 'd'],bk);
            jj1 = jj2+1;
            jj2 = jj2+npts(kk);
            mm1 = mm2+1;
            mm2 = mm2 + nverts(kk);
            XXk = XX(:,jj1:jj2);
            vertsk = vertices(:,mm1:mm2);
            xxk = xx(:, jj1:jj2)+1;
            verts2dk = verts2d(:, mm1:mm2)+1;
            subface = reshape(repmat(cface(:),1,npts(kk))+repmat(8*(0:npts(kk)-1),24,1), 4, []);
            % draw 2D projection in all views
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
            text(tpx,tpy,['\it\fontname{Arial}\fontsize{9}\color{white}Points - ' ...
                '\color{yellow}',framenb],'HorizontalAlignment','right');
            for pp = 1:n_cam,
                if active_imgviews(pp,bk),
                    ii = (pp-1)*2;
                    x_kk = xxk(ii+1:ii+2, :);
                    plot(x_kk(1,:),x_kk(2,:),'g.');
                    if draw_faces,
                        x_kk = verts2dk(ii+1:ii+2, :);
                        patch('Faces', subface', 'Vertices', x_kk', 'FaceColor', 'b', 'EdgeColor','r','FaceAlpha',0.05);
                    end;
                end;
            end;
            % set axis position in window, position = [left, bottom, width, height]
            set(gca,'position',[0 0 1 1]);
            set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
            drawnow;
            save_name = [imgdir '/map_' framenb '.jpg'];
            print(2,'-djpeg',['-r' num2str(resolution)],save_name);
            
            % plot 3D points on border and perspective image in 3D space
            figure(3); hold off;
            plot3(XXk(1,:),XXk(3,:),-XXk(2,:),'g.');
            hold on;
            if pim,
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
                        XXk = hc+img_position*vect/norm(vect);
                        % position of pespective image corners (o, y, xy, x)
                        XXk = XXk(:,ones(1,4))+(Rkk(:,1)*[-1,-1,1,1]+Rkk(:,2)*[-1,1,1,-1])*maxbound*0.6;
                        % project the four corner on image
                        x_kk = project_points_mirror2(XXk,om,T,hand,fc,cc,zeros(5,1),alpha_c);
                        Ikk = homography_image(frame_kk, x_kk+1, [], 0);   % image background is black
                        % transform to matlab 3D axes for surface function (o, y, x, xy)
                        XXk = XXk(:, [1,2,4,3]);
                        xc(:) = XXk(1,:);
                        yc(:) = XXk(2,:);
                        zc(:) = XXk(3,:);
                        surface(xc, zc, -yc, Ikk,'FaceColor','texturemap','EdgeColor','none');
                    end;
                end;
            end;
            % set 3D plot
            figure(3);
            if draw_faces,
                patch('Faces', subface', 'Vertices', [vertsk(1,:); vertsk(3,:); -vertsk(2,:)]', 'FaceColor', 'b','EdgeColor','r','FaceAlpha',0.1);
                if Rd,
                    set(3, 'renderer', 'zbuffer'); cameratoolbar('ResetCameraAndSceneLight');
                    cameratoolbar('Togglescenelight');
                end;
            end;
            axis equal tight vis3d off;
            view(az,el);
            set(3,'color',[0 0 0], 'PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
            drawnow;
            save_name = [imgdir '/hull_' framenb '.jpg'];
            print(3,'-djpeg',['-r' num2str(resolution)],save_name);
            saveas(3, [imgdir '/hull_' framenb],'fig');
        end;
        
    else
        
        xx = NaN(2*n_cam, sum_npts);
        verts2d = NaN(2*n_cam, sum_npts*8);
        for pp = 1:n_cam,
            fc = fc_mat(:,pp);
            cc = cc_mat(:,pp);
            alpha_c = alpha_vec(pp);
            om = Omcc(:,pp);
            T = Tcc(:,pp);
            hand = handcc(pp);
            x_kk = project_points_mirror2([XX,vertices],om,T,hand,fc,cc,zeros(5,1),alpha_c);
            ii = (pp-1)*2;
            xx(ii+1:ii+2, :) = x_kk(:,1:sum_npts);
            verts2d(ii+1:ii+2, :) = x_kk(:,sum_npts+1:end);
        end;
        
        jj2 = 0; mm2 = 0;
        for kk = ind_active,
            bk = base+kk;
            framenb = sprintf(['%0' ndigit 'd'],bk);
            jj1 = jj2+1;
            jj2 = jj2+npts(kk);
            mm1 = mm2+1;
            mm2 = mm2 + nverts(kk);
            XXk = XX(:,jj1:jj2);
            vertsk = vertices(:,mm1:mm2);
            xxk = xx(:, jj1:jj2)+1;
            verts2dk = verts2d(:, mm1:mm2)+1;
            subface = reshape(repmat(cface(:),1,npts(kk))+repmat(8*(0:npts(kk)-1),24,1), 4, []);
            % plot visual hull and perspective image in 3D space
            figure(3); hold off;
            plot3(XXk(1,:),XXk(3,:),-XXk(2,:),'g.');
            hold on;
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
                    x_kk = xxk(ii+1:ii+2, :);
                    plot(x_kk(1,:),x_kk(2,:),'g.');
                    if draw_faces,
                        x_kk = verts2dk(ii+1:ii+2, :);
                        patch('Faces', subface', 'Vertices', x_kk', 'FaceColor', 'b', 'EdgeColor','r','FaceAlpha',0.05);
                    end;
                    % set axis position in window, position = [left, bottom, width, height]
                    set(gca,'position',[0 0 1 1]);
                    set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                    drawnow;
                    save_name = [imgdir_cell{pp} '/map_cam' num2str(pp) '_' framenb '.jpg'];
                    print(2,'-djpeg',['-r' num2str(resolution)],save_name);
                    % draw 3D
                    if pim,
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
                        XXk = hc+img_position*vect/norm(vect);
                        % position of pespective image corners (o, y, xy, x)
                        XXk = XXk(:,ones(1,4))+(Rkk(:,1)*[-1,-1,1,1]+Rkk(:,2)*[-1,1,1,-1])*maxbound*0.6;
                        % project the four corner on image
                        x_kk = project_points_mirror2(XXk,om,T,hand,fc,cc,zeros(5,1),alpha_c);
                        Ikk = homography_image(frame_kk, x_kk+1, [], 0);   % image background is black
                        % transform to matlab 3D axes for surface function (o, y, x, xy)
                        XXk = XXk(:, [1,2,4,3]);
                        xc(:) = XXk(1,:);
                        yc(:) = XXk(2,:);
                        zc(:) = XXk(3,:);
                        figure(3);
                        surface(xc, zc, -yc, Ikk,'FaceColor','texturemap','EdgeColor','none');
                    end;
                end;
            end;
            % set 3D plot
            figure(3);
            if draw_faces,
                patch('Faces', subface', 'Vertices', [vertsk(1,:); vertsk(3,:); -vertsk(2,:)]', 'FaceColor', 'b','EdgeColor','r','FaceAlpha',0.1);
                if Rd,
                    set(3, 'renderer', 'zbuffer'); cameratoolbar('ResetCameraAndSceneLight');
                    cameratoolbar('Togglescenelight');
                end;
            end;
            axis equal tight vis3d off;
            view(az,el);
            nx = imsize(1,1);
            ny = imsize(2,1);
            resolution = round(ny/height);
            set(3,'color',[0 0 0], 'PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
            drawnow;
            % set the saving directory to the first camera view folder
            imgdir = imgdir_cell{1};
            save_name = [imgdir '/hull_' framenb '.jpg'];
            print(3,'-djpeg',['-r' num2str(resolution)],save_name);
            saveas(3, [imgdir '/hull_' framenb],'fig');
        end;
    end;
end;
ind_active = find(active_images);
