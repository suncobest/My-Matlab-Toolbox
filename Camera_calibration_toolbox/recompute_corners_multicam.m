% Re-extract image pattern corners after calibration
if ~exist('imstrnum','var')||~exist('fc_mat','var'),
    fprintf(1,'No calibration data available.\n');
    return;
end;

ind_active_views = find(active_imgviews(:)');
if isempty(ind_active_views),
    fprintf(1,'No image is available!\n');
    return;
end;

flag = 0;
if exist('y_cell','var'),
    for kk = ind_active_views,
        x_kk = x_cell{kk};
        y_kk = y_cell{kk};
        if isempty(y_kk) || isnan(any(y_kk(:))) || isempty(x_kk) || isnan(any(x_kk(:))),
            flag = 1;
            fprintf(1,'\nWARNING: Data of x_cell{%d} or y_cell{%d} not available! Check it out!\n',kk,kk);
            break;
        end;
    end;
end;

if ~exist('imsize','var'),
    fprintf(1,'WARNING: No image size (imsize) available!''\n');
    flag = 1;
end;

if flag,
    fprintf(1,['Need to calibrate once before recomputing image corners.\n' ...
        'You can load ''Multicam_Calib_Results.mat'' if the file is available ...\n']);
    return;
end;

if ~exist('win_size', 'var'),
    win_size = zeros(2,n_cam);
    for pp=1:n_cam,
        fprintf(1,'\nWindow size for corner finder (width and height) of view %d:\n', pp);
        wintx = input(['wintx = ([] = ' num2str(wintx_default) ')']);
        if isempty(wintx), wintx = wintx_default; end;
        wintx = round(wintx);
        winty = input(['winty = ([] = ' num2str(winty_default) ')']);
        if isempty(winty), winty = winty_default; end;
        winty = round(winty);
        win_size(:, pp) = [wintx; winty];
        fprintf(1,'Window size of view %d = %dx%d\n', pp,2*wintx+1,2*winty+1);
    end;
end;

fprintf(1,'\nRe-extraction of the grid corners on the images (after first calibration):\n');
% Color code for each image:
palette = 'brgkcm';
% set waiting time to refresh images
waitDt = 0.3;
Txt_shift = max(imsize,1)/40;         % text offset
ind_active = find(active_images);
flag = length(ind_active)*n_cam>30;
no_grid = 0;
if ~exist('n_sq_mat','var'),
    no_grid = 1;
end;

for kk = ind_active,
    rgbi = palette(rem(kk-1,6)+1);
    active_view = active_imgviews(:,kk);
    for pp =1:n_cam,
        if active_view(pp),
            format_image = imformat{pp};
            ima_name = [imbase{pp} imstrnum{pp}{kk} '.' format_image];
            if exist(ima_name,'file')==2,
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
                % display the corners on image
                figure(6);
                image(I); axis image;
                colormap(gray(256));
                title(['(Camera ' num2str(pp) ' Image ' num2str(kk) ') - New corner points (+) and old points (o)']);
                hold on;

                kth = (kk-1)*n_cam+pp;
                x_kk = x_cell{kth};     % old corner points
                y_kk = y_cell{kth};     % projected points
                % recompute corners given projection
                y_kk = cornerfinder(y_kk+1,I,win_size(1,pp),win_size(2,pp));
                figure(6);
                plot(x_kk(1,:)+1,x_kk(2,:)+1,[rgbi, 'o'],y_kk(1,:),y_kk(2,:),'r+');
                % refresh corner points
                x_kk = y_kk-1;
                x_cell{kth} = x_kk;
                % plot axis on the figure:
                if ~no_grid,
                    N_kk = size(x_kk,2);
                    n_sq_x = n_sq_mat(1,kth);
                    n_sq_y = n_sq_mat(2,kth);
                    if (N_kk ~= ((n_sq_x+1)*(n_sq_y+1))),
                        no_grid = 1;
                    end;
                    n_sq_x = n_sq_mat(1,kth);
                    n_sq_y = n_sq_mat(2,kth);

                    Nx = n_sq_x+1;
                    Ny = n_sq_y+1;

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

                    plot([xo;xX]+1,[yo;yX]+1,'g-','linewidth',2);
                    plot([xo;xY]+1,[yo;yY]+1,'g-','linewidth',2);
                    delta = Txt_shift(pp);
                    text(xXc + delta * uX(1) +1, yXc + delta * uX(2)+1,'X','color','g','Fontsize',14,'HorizontalAlignment','center');
                    text(xYc + delta * uY(1)+1, yYc + delta * uY(2)+1,'Y','color','g','Fontsize',14,'HorizontalAlignment','center');
                    text(xo + delta * uO(1) +1, yo + delta * uO(2)+1,'O','color','g','Fontsize',14,'HorizontalAlignment','center');
                end;
                % set figure 6
                set(6,'color',[1 1 1],'Name',num2str(kk),'NumberTitle','off');
                zoom on;
                drawnow;
                hold off;
                if flag,
                    pause(waitDt);
                else
                    fprintf(1,['\nCheck the recomputed new corners in (camera %d, image %d), ' ...
                        'press any key to continue...\n'], pp, kk);
                    pause;
                end;
            else
                fprintf(1,'NOTE: Image %s not found!\n',ima_name);
            end;
        end;
    end;
end;

% Recompute the error:
y_cell = cell(1, n_view);  % Reprojected points
ex_cell = y_cell;             % Reprojected error
err_cam = zeros(2,n_cam);
ex = [];
if exist('Qw_mat','var'),
    for pp = 1:n_cam,
        fc = fc_mat(:,pp);
        cc = cc_mat(:,pp);
        alpha_c = alpha_vec(pp);
        kc = kc_mat(:,pp);
        handkk = hand_list(pp);
        % Calibration matrix:
        ind_active_views = find(active_imgviews(pp,:));
        for kk = ind_active_views,
            kth = (kk-1)*n_cam+pp;
            Qwkk = Qw_mat(:, kth);
            Twkk = Tw_mat(:, kth);
            X_kk = X_cell{kth};
            x_kk = x_cell{kth};
            y_kk = project_points_mirror(X_kk,Qwkk,Twkk,handkk,fc,cc,kc,alpha_c);
            ex_kk = x_kk-y_kk;
            ex = [ex, ex_kk];
            y_cell{kth} = y_kk;
            ex_cell{kth} = ex_kk;
        end;
        err_cam(:,pp) = std(cell2mat(ex_cell(pp:n_cam:end)),0,2);
    end;
else
    for pp = 1:n_cam,
        fc = fc_mat(:,pp);
        cc = cc_mat(:,pp);
        alpha_c = alpha_vec(pp);
        kc = kc_mat(:,pp);
        handkk = hand_list(pp);
        % Calibration matrix:
        ind_active_views = find(active_imgviews(pp,:));
        for kk = ind_active_views,
            kth = (kk-1)*n_cam+pp;
            omwkk = Omw_mat(:, kth);
            Twkk = Tw_mat(:, kth);
            X_kk = X_cell{kth};
            x_kk = x_cell{kth};
            y_kk = project_points_mirror2(X_kk,omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
            ex_kk = x_kk-y_kk;
            ex = [ex, ex_kk];
            y_cell{kth} = y_kk;
            ex_cell{kth} = ex_kk;
        end;
        err_cam(:,pp) = std(cell2mat(ex_cell(pp:n_cam:end)),0,2);
    end;
end;
err_std = std(ex,0,2);
fprintf(1,'\ndone\n');
