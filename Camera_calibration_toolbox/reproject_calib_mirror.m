%%%%%%%%%%%%%%%%%%%% REPROJECT ON THE IMAGES %%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('fc','var')||~exist('y_cell','var'),
    fprintf(1,'No calibration data available.\n');
    return;
end;

if ~exist('nx','var')||~exist('ny','var'),
    fprintf(1,'WARNING: No image size (nx,ny) available. Setting nx=640 and ny=480\n');
    nx = 640;
    ny = 480;
end;
delta = max(nx, ny)/40;        % text offset

% Color code for each image:
palette = 'brgkcm';
% set waiting time to refresh images
waitDt = 0.3;

% Reproject the patterns on the images, and compute the pixel errors:
% Reload the images if available
ind_active_views = find( active_imgviews(1,:));
if isempty(ind_active_views),
    fprintf(1,'No image of view 1 is available! Please check your data!\n');
    return;
end;

% draw the reprojection of points
flag = length(ind_active_views)>20;
no_grid = 0;
if ~exist('n_sq_mat','var'),
    no_grid = 1;
end;

for kk = ind_active_views,
    string_num = strnum_cell{kk};                           % strnum_cell为按顺序排列的图片序号字符cell
    ima_name = [calib_name  string_num '.' format_image];
    if exist(ima_name,'file')~=2,
        I = 255*ones(ny,nx);
    else
        fprintf(1,'Loading image %s...\n',ima_name);
        if format_image(1) == 'p',
            if format_image(2) == 'p',
                I = double(loadppm(ima_name));
            elseif format_image(2) == 'g',
                I = double(loadpgm(ima_name));
            else
                I = double(imread(ima_name));
            end;
        else
            if format_image(1) == 'r',
                I = readras(ima_name);
            else
                I = double(imread(ima_name));
            end;
        end;

        if size(I,3)==3,
            I = 0.299 * I(:,:,1) + 0.5870 * I(:,:,2) + 0.114 * I(:,:,3);
        end;
    end;
    figure(6);
    image(I); axis image;
    colormap(gray(256));
    title(['Image ' num2str(kk) ' - Image points (+) and reprojected grid points (o)']);
    hold on;
    rgbi = palette(rem(kk-1,6)+1);
    active_view = active_imgviews(:,kk);
    for pp = 1:n_cam,
        if active_view(pp),
            kth = (kk-1)*n_cam+pp;
            x_kk = x_cell{kth};
            y_kk = y_cell{kth};
            if isempty(x_kk) || isempty(y_kk),
                 fprintf(1,'\nWARNING: Data of view %d of frame %d not available!\n',pp,kk);
                continue;
            end;
            figure(6);
            plot(x_kk(1,:)+1,x_kk(2,:)+1,'r+',y_kk(1,:)+1,y_kk(2,:)+1,[rgbi, 'o']);
            if ~no_grid,
                N_kk = size(x_kk,2);
                Nx = n_sq_mat(1,kth)+1;
                Ny = n_sq_mat(2,kth)+1;
                assert(N_kk == Nx*Ny,'Number of corner points do not match with square size!');
                % plot more things on the figure (to help the user):
                ind_ori = (Ny - 1) * Nx + 1;
                ind_X = Nx*Ny;
                ind_Y = 1;
                ind_XY = Nx;

                xo = x_kk(1,ind_ori);
                yo = x_kk(2,ind_ori);

                xX = x_kk(1,ind_X);
                yX = x_kk(2,ind_X);

                xY = x_kk(1,ind_Y);
                yY = x_kk(2,ind_Y);

                xXY = x_kk(1,ind_XY);
                yXY = x_kk(2,ind_XY);

                % center of grid
                uu = cross(cross([xo;yo;1],[xXY;yXY;1]),cross([xX;yX;1],[xY;yY;1]));
                xc = uu(1)/uu(3);
                yc = uu(2)/uu(3);

                % the joint point of X axis and the joint line of Y-directional
                % infinite point and the center (过中心沿-y方向与x轴的交点)
                bbb = cross(cross([xo;yo;1],[xY;yY;1]),cross([xX;yX;1],[xXY;yXY;1]));
                uu = cross(cross([xo;yo;1],[xX;yX;1]),cross([xc;yc;1],bbb));
                xXc = uu(1)/uu(3);
                yXc = uu(2)/uu(3);

                % the joint point of Y axis and the joint line of X-directional
                % infinite point and the center (过中心沿-x方向与y轴的交点)
                bbb = cross(cross([xo;yo;1],[xX;yX;1]),cross([xY;yY;1],[xXY;yXY;1]));
                uu = cross(cross([xo;yo;1],[xY;yY;1]),cross([xc;yc;1],bbb));
                xYc = uu(1)/uu(3);
                yYc = uu(2)/uu(3);

                uX = [xXc - xc;yXc - yc];
                uY = [xYc - xc;yYc - yc];
                uO = [xo - xc;yo - yc];

                uX = uX / norm(uX);
                uY = uY / norm(uY);
                uO = uO / norm(uO);

                figure(6);
                plot([xo;xX]+1,[yo;yX]+1,'g-','linewidth',2);
                plot([xo;xY]+1,[yo;yY]+1,'g-','linewidth',2);
                text(xXc + delta * uX(1) +1, yXc + delta * uX(2)+1,'X','color','g','Fontsize',14,'HorizontalAlignment','center');
                text(xYc + delta * uY(1)+1, yYc + delta * uY(2)+1,'Y','color','g','Fontsize',14,'HorizontalAlignment','center');
                text(xo + delta * uO(1) +1, yo + delta * uO(2)+1,'O','color','g','Fontsize',14,'HorizontalAlignment','center');
            end;
        end;
    end;
    % set figure 6
    set(6,'color',[1 1 1],'Name',num2str(kk),'NumberTitle','off');
    zoom on;
    drawnow;
    hold off;
    if flag,
        pause(waitDt);
    else
        fprintf(1,['\nCheck the reprojected corners on image %d, ' ...
            'press any key to continue...\n'], kk);
        pause;
    end;
end;
fprintf(1,'Pixel error:    err = [%3.5f   %3.5f] (all active image views)\n\n',err_std);
