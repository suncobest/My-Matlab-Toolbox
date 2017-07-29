% Polygon projection: use wing center to set wing's position, not wing root

% column number limit of each member
Nb = [1, nptb+3; nptb+4, nptb+nptw+6; nptb+nptw+7, nptb+2*nptw+9];
% faces of every member
XYZs = {XYZb, XYZw, XYZw.*([1;-1;1]*ones(1,nptw))};
faces = {faceb, facew, facew};
npts = [nptb,nptw,nptw];
mainAx = [axes_mean(1,1), axes_mean(2,2), axes_mean(2,3)];

if calib_mode,
    frame_kk = imread([imgdir '/' imgbase strnum_frame{1} '.' imgfmt]);
    SW = size(frame_kk,3)==3;
    resolution = round(ny/height);   % dpi
    for kk=ind_active,
        XX = NaN(3,Nb(end));
        % plot 3D axes and polygon
        figure(4); hold off;
        for ii = 1:nparts+1,
            ax = axes_mean(:,ii);
            ct = member_center(:,ii,kk);
            vc = rodrigues(member_omega(:,ii,kk));
            ctt = ct(:,ones(1,2));
            xyz = ctt+vc(:,1:2)*diag(ax(1:2))*nlen_vec;
            % resize, reorient, and relocate polygon
            if norm_sw,
                XYZe = vc*mainAx(ii)*XYZs{ii}+ct(:,ones(1,npts(ii)));
            else
                XYZe = vc*diag(ax)*XYZs{ii}+ct(:,ones(1,npts(ii)));
            end;
            XX(:, Nb(ii,1) : Nb(ii,2)) = [ct,xyz,XYZe];
            for i=1:2,
                plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',w3);
                hold on;
            end;
            patch('Faces', faces{ii}, 'Vertices', [XYZe(1,:); XYZe(3,:); -XYZe(2,:)]','FaceColor',palette(ii,:),...
                'EdgeColor',palette2(ii,:),'FaceAlpha',al);
        end;

        % superimpose 2D axes and polygon on the original image (figure 2)
        framenb = strnum_frame{kk};
        frame_kk = imread([imgdir '/' imgbase framenb '.' imgfmt]);
        if SW,
            frame_kk = 0.299 * frame_kk(:,:,1) + 0.587 * frame_kk(:,:,2) + 0.114 * frame_kk(:,:,3);
        end;
        figure(2); hold off;
        image(frame_kk);
        colormap(map);      % colormap for gray scale image
        hold on;
        % generate synthetic silhouette (figure 3)
        frame_kk = false(ny,nx);
        figure(3); hold off;
        image(frame_kk);
        colormap(gray(2));      % colormap for binary image
        hold on;
        for pp = 1:n_cam,
            if active_imgviews(pp,kk),
                om = Omcc(:,pp);
                T = Tcc(:,pp);
                hand = handcc(pp);
                xx = project_points_mirror2(XX,om,T,hand,fc,cc,zeros(5,1),alpha_c);
                for ii =1:nparts+1,
                    ctt = repmat(xx(:, Nb(ii,1)),1,2);
                    xyz = xx(:, Nb(ii,1)+1 : Nb(ii,1)+2);
                    XYZe = xx(:, Nb(ii,1)+3 : Nb(ii,2));
                    figure(2);
                    for i=1:2,
                        plot([ctt(1,i);xyz(1,i)], [ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',w1);
                    end;
                    patch('Faces', faces{ii}, 'Vertices', XYZe'+1, 'FaceColor',palette(ii,:),'EdgeColor',palette2(ii,:),...
                        'LineWidth',w2,'FaceAlpha',al);
                    figure(3);
                    patch('Faces', faces{ii}, 'Vertices', XYZe'+1, 'FaceColor','w','EdgeColor','none');
                end;
            end;
        end;
        figure(2);
        set(gca,'position',[0 0 1 1]); axis image off;
        set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
        drawnow;
        save_name = [imgdir '/projection_' imgbase framenb '.jpg'];
        print(2,'-djpeg',['-r' num2str(resolution)],save_name);
        figure(3);
        set(gca,'position',[0 0 1 1]);
        set(3,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
        drawnow;
        temp = [imgdir '/silhouette_' imgbase framenb '.png'];
        print(3,'-dpng',['-r' num2str(resolution)], temp);
        frame_kk = imread(temp);
        frame_kk = im2bw(frame_kk(:,:,1),0.5);
        imwrite(frame_kk,temp,'png');
        
        % plot perspective image in 3D space
        figure(4);
        frame_kk = imread(save_name);
        for pp = 1:n_cam,
            if active_imgviews(pp,kk),
                xx = xim((pp-1)*2+1 : pp*2, :);
                Ikk = homography_image(frame_kk, xx+1);   % image background is white
                XX = Xim((pp-1)*3+1 : pp*3, :);
                xc(:) = XX(1,:);
                yc(:) = XX(2,:);
                zc(:) = XX(3,:);
                surface(xc, zc, -yc, Ikk,'FaceColor','texturemap','EdgeColor','none');
            end;
        end;
        % set 3D plot
        axis equal tight vis3d off;
        set(4, 'renderer', 'zbuffer','color',[1,1,1]*0.7);
        set(4, 'PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
        cameratoolbar('ResetCameraAndSceneLight');
        cameratoolbar('Togglescenelight');
        view(az,el);
        if kk==ind_active(1),
            re_edit = 1;
            while re_edit,
                temp = input('Set looking direction for 3D view: ([azimuth, elevation]) ');
                if length(temp)~=2,
                    fprintf(1,'Unexpected input! Please enter again!\n');
                    continue;
                end;
                az = temp(1);
                el = temp(2);
                figure(4);
                view(az,el);
                re_edit = input('Need to reset looking direction or not? ([]=yes, other=no) ','s');
                re_edit = isempty(re_edit);
            end;
        end;
        if set_range,
            axis(axes_range);
        end;
        drawnow;
        save_name = [imgdir '/3D_projection_' framenb '.jpg'];
        print(4,'-djpeg',['-r' num2str(resolution)],save_name);
        saveas(4, [imgdir '/3D_projection_' framenb],'fig');
    end;
    
else
    
    for kk=ind_active,
        XX = NaN(3,Nb(end));
        % plot 3D axes and polygon
        figure(4); hold off;
        for ii = 1:nparts+1,
            ax = axes_mean(:,ii);
            ct = member_center(:,ii,kk);
            vc = rodrigues(member_omega(:,ii,kk));
            ctt = ct(:,ones(1,2));
            xyz = ctt+vc(:,1:2)*diag(ax(1:2))*nlen_vec;
            % resize, reorient, and relocate polygon
            if norm_sw,
                XYZe = vc*mainAx(ii)*XYZs{ii}+ct(:,ones(1,npts(ii)));
            else
                XYZe = vc*diag(ax)*XYZs{ii}+ct(:,ones(1,npts(ii)));
            end;
            XX(:, Nb(ii,1) : Nb(ii,2)) = [ct,xyz,XYZe];
            for i=1:2,
                plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',w3);
                hold on;
            end;
            hold on;
            patch('Faces', faces{ii}, 'Vertices', [XYZe(1,:); XYZe(3,:); -XYZe(2,:)]','FaceColor',palette(ii,:),...
                'EdgeColor',palette2(ii,:),'FaceAlpha',al);
        end;
  
        for pp = 1:n_cam,
            fc = fc_mat(:,pp);
            cc = cc_mat(:,pp);
            alpha_c = alpha_vec(pp);
            om = Omcc(:,pp);
            T = Tcc(:,pp);
            hand = handcc(pp);
            xx = project_points_mirror2(XX,om,T,hand,fc,cc,zeros(5,1),alpha_c);
            if active_imgviews(pp,kk),
                imgdir = imgdir_cell{pp};
                imgbase = imgbase_cell{pp};
                framenb = strnum_frameNcam{pp}{kk};
                % superimpose 2D axes and ellipsoid on the original image (figure 2)
                frame_kk = imread([imgdir '/' imgbase framenb '.' imgfmt_cell{pp}]);
                if size(frame_kk,3)==3,
                    frame_kk = 0.299 * frame_kk(:,:,1) + 0.587 * frame_kk(:,:,2) + 0.114 * frame_kk(:,:,3);
                end;
                figure(2); hold off;
                image(frame_kk);
                colormap(map);      % colormap for gray scale image
                hold on;
                % generate synthetic silhouette (figure 3)
                nx = imsize(1,pp);
                ny = imsize(2,pp);
                resolution = round(ny/height);   % dpi
                frame_kk = false(ny,nx);
                figure(3); hold off;
                image(frame_kk);
                colormap(gray(2));      % colormap for binary image
                hold on;
                for ii = 1:nparts+1,
                    ctt = repmat(xx(:, Nb(ii,1)),1,2);
                    xyz = xx(:, Nb(ii,1)+1 : Nb(ii,1)+2);
                    XYZe = xx(:, Nb(ii,1)+3 : Nb(ii,2));
                    figure(2);
                    for i=1:2,
                        plot([ctt(1,i);xyz(1,i)], [ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',w1);
                    end;
                    patch('Faces', faces{ii}, 'Vertices', XYZe'+1, 'FaceColor',palette(ii,:),'EdgeColor',palette2(ii,:),...
                        'LineWidth',w2,'FaceAlpha',al);
                    figure(3);
                    patch('Faces', faces{ii}, 'Vertices', XYZe'+1, 'FaceColor','w','EdgeColor','none');
                end;
                figure(2);
                set(gca,'position',[0 0 1 1]); axis image off;
                set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                drawnow;
                save_name = [imgdir '/projection_' imgbase framenb '.jpg'];
                print(2,'-djpeg',['-r' num2str(resolution)],save_name);
                figure(3);
                set(gca,'position',[0 0 1 1]);
                set(3,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                drawnow;
                temp = [imgdir '/silhouette_' imgbase framenb '.png'];
                print(3,'-dpng',['-r' num2str(resolution)], temp);
                frame_kk = imread(temp);
                frame_kk = im2bw(frame_kk(:,:,1),0.5);
                imwrite(frame_kk,temp,'png');
                % draw 3D: plot perspective image in 3D space
                frame_kk = imread(save_name);
                xx = xim((pp-1)*2+1 : pp*2, :);
                Ikk = homography_image(frame_kk, xx+1);   % image background is white
                Xk = Xim((pp-1)*3+1 : pp*3, :);
                xc(:) = Xk(1,:);
                yc(:) = Xk(2,:);
                zc(:) = Xk(3,:);
                figure(4);
                surface(xc, zc, -yc, Ikk,'FaceColor','texturemap','EdgeColor','none');
            end;
        end;
        
        % set 3D plot
        figure(4);
        axis equal tight vis3d off;
        set(4, 'renderer', 'zbuffer','color',[1,1,1]*0.7);
        set(4, 'PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
        cameratoolbar('ResetCameraAndSceneLight');
        cameratoolbar('Togglescenelight');
        view(az,el);
        if kk==ind_active(1),
            re_edit = 1;
            while re_edit,
                temp = input('Set looking direction for 3D view: ([azimuth, elevation]) ');
                if length(temp)~=2,
                    fprintf(1,'Unexpected input! Please enter again!\n');
                    continue;
                end;
                az = temp(1);
                el = temp(2);
                figure(4);
                view(az,el);
                re_edit = input('Need to reset looking direction or not? ([]=yes, other=no) ','s');
                re_edit = isempty(re_edit);
            end;
        end;
        if set_range,
            axis(axes_range);
        end;
        drawnow;
        imgdir = imgdir_cell{1};
        framenb = strnum_frameNcam{1}{kk};
        save_name = [imgdir '/3D_projection_' framenb '.jpg'];
        print(4,'-djpeg',['-r' num2str(resolution)],save_name);
        saveas(4, [imgdir '/3D_projection_' framenb],'fig');
    end;
end;
