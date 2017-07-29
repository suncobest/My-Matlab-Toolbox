if ~exist('fc_mat','var')||~exist('y_cell','var'),
    fprintf(1,'No calibration data available.\n');
    return;
end;
ind_active_views = find(active_imgviews(1,:));
if isempty(ind_active_views),
    fprintf(1,'No image views are activated! Please check your data!\n');
    return;
end;

% Color code for each image:
palette = 'brgkcm';

no_grid = 0;
if ~exist('n_sq_mat','var'),
    no_grid = 1;
end;
% draw the distribution of error
figure(5);
hold on;
for kk = ind_active_views,
    rgbi = palette(rem(kk-1,6)+1);
    active_view = active_imgviews(:,kk);
    for pp = 1:n_cam,
        if active_view(pp),
            kth = (kk-1)*n_cam+pp;
            ex_kk = ex_cell{kth};
            if isempty(ex_kk),
                fprintf(1,'\nWarning: the error of (camera %d, image %d) not available! Please check your data!\n',pp,kk);
                continue;
            end;
            plot(ex_kk(1,:),ex_kk(2,:),[rgbi,'+']);
        end;
    end;
end;
set(5,'color',[1 1 1],'Name','error','NumberTitle','off');
title('Reprojection error (in pixel) - To exit: right button');
xlabel('x');
ylabel('y');
axis equal;
drawnow;
hold off;

fprintf(1,'Pixel error:      err = [%3.5f   %3.5f] (all active images)\n\n',err_std);

ex = cell2mat(ex_cell);
lmb = 1;      % ginput "left mouse button": 1=LMB, 2=MMB, 3=RMB
[xp,yp] = ginput(1);
ind_active_views = find(active_imgviews(:)');
while lmb==1,
    ddd = (ex(1,:)-xp).^2 + (ex(2,:)-yp).^2;
    [~,indmin] = min(ddd);
    ex_mx = ex(1,indmin);
    ex_my = ex(2,indmin);
    for kth = ind_active_views,
        ex_kk = ex_cell{kth};
        sol_kk = find((ex_kk(1,:) == ex_mx)&(ex_kk(2,:) == ex_my));
        if ~isempty(sol_kk),
            break;
        end;
    end;
    x_kk = x_cell{kth};
    xpt = x_kk(:,sol_kk);
    kk = floor((kth-1)/n_cam)+1;
    pp = kth-(kk-1)*n_cam;
    if ~no_grid,
        n_sq_x = n_sq_mat(1,kth);
        n_sq_y = n_sq_mat(2,kth);
        Nx = n_sq_x+1;
        Ny = n_sq_y+1;
        
        y1 = floor((sol_kk-1)/Nx);
        x1 = sol_kk  - Nx*y1;   % the 2D index (i,j)
        y1 = Ny - y1;
        
        fprintf(1,'\nSelected point from: (camera %d, image %d)\n',pp,kk);
        fprintf(1,'Selected point index: %d\n',sol_kk);
        fprintf(1,'Pattern coordinates (in units of (dX,dY)): (X,Y)=(%d,%d)\n',[x1-1 y1-1]);
        fprintf(1,'Image coordinates (in pixel): (%3.2f,%3.2f)\n',xpt);
        fprintf(1,'Pixel error = (%3.5f,%3.5f)\n',[ex_mx ex_my]);
    else
        fprintf(1,'\nSelected point from: (camera %d, image %d)\n',pp,kk);
        fprintf(1,'Selected point index: %d\n',sol_kk);
        fprintf(1,'Image coordinates (in pixel): (%3.2f,%3.2f)\n',xpt);
        fprintf(1,'Pixel error = (%3.5f,%3.5f)\n',[ex_mx ex_my]);
    end;
    [xp,yp,lmb] = ginput(1);
end;

disp('done');
