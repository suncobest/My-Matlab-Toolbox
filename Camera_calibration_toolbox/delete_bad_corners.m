% delete one row or one column of corner points
if ~exist('x_cell','var'),
    fprintf(1,'\nThere is no corner points data to process!\n');
    return;
end;
if ~exist('y_cell','var'),
    fprintf(1,'\nThere is no reprojected corner points data! Please calibrate first!\n');
    return;
end;
fprintf(1,'\nThis script could delete bad rows or columns of a specific corner matrix...\n');
if ~exist('pp','var') || ~exist('kk','var'),
    flag = 1;
else
    flag = ~active_imgviews(pp,kk);
end;
while flag,
    fprintf(1,'You have to specify active corner points:\n');
    pp = input(['Specify the camera view of corner data to process: ([' num2str(1:n_cam) '])']);
    assert(length(pp)==1 && any(pp==(1:n_cam)),'Unexpected input! Please choose one view at a time!');
    kk = input(['Specify the image number of corner data to process: (1 ~' num2str(n_ima) ')']);
    assert(length(kk)==1 && any(kk==(1:n_ima)),'Unexpected input! Please indicate one image!');
    flag = ~active_imgviews(pp,kk);
end;

calib_name = imbase{pp};
format_image = imformat{pp};
strnum_cell = imstrnum{pp};
string_num = strnum_cell{kk};      % strnum_cell contain sorted string numbers of images
ima_name = [calib_name  string_num '.' format_image];
if exist(ima_name,'file')==2,
    fprintf(1,'Processing with (camera %d, image %d)...\n',pp,kk);
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
end
kth = (kk-1)*n_cam+pp;
x_kk = x_cell{kth};
y_kk = y_cell{kth};
X_kk=X_cell{kth};
% plot axes and corner points on image
figure(6),image(I);
colormap(map);
hold on;
N_kk = size(x_kk,2);
Nx = n_sq_mat(1,kth)+1;
Ny = n_sq_mat(2,kth)+1;
assert(N_kk == Nx*Ny,'Number of corner points do not match with square size!');

ind_ori = (Ny - 1) * Nx + 1;
ind_X = N_kk;
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
nx = imsize(1,pp);
ny = imsize(2,pp);
delta = max(nx, ny)/40;        % text offset
plot([xo;xX]+1,[yo;yX]+1,'g-','linewidth',1);
plot([xo;xY]+1,[yo;yY]+1,'g-','linewidth',1);
plot(y_kk(1,:)+1,y_kk(2,:)+1,'go');
plot(x_kk(1,:)+1,x_kk(2,:)+1,'r+','markersize',10);
text(xXc + delta * uX(1) +1, yXc + delta * uX(2)+1,'X','color','g','Fontsize',14,'HorizontalAlignment','center');
text(xYc + delta * uY(1)+1, yYc + delta * uY(2)+1,'Y','color','g','Fontsize',14,'HorizontalAlignment','center');
text(xo + delta * uO(1) +1, yo + delta * uO(2)+1,'O','color','g','Fontsize',14,'HorizontalAlignment','center');
title(['Check the extracted corners of (camera ' num2str(pp) ', image ' num2str(kk) '):']);
hold off;

flag = input('Is this image view good enough to keep active? ([]=yes, other=no) ','s');
if ~isempty(flag),
    active_imgviews(pp,kk) = 0;
    fprintf(1,'(camera %d, image %d) is now set inactive.\n',pp,kk);
    return;
end;

flag = input('Need to delete bad rows or columns? ([]=no, other=yes) ','s');
if isempty(flag),
    fprintf(1,'Nothing is changed by now!\n');
    return;
end;

%% corner matrix (delete from xy, cannot delete o)
%  o(3)________________________y(1)
%     |           |           |
%     |           |           |
%     |--------------------------row
%     |           |           |
%     |___________|___________|
%    x(4)         |            xy(2)
%             column


flag = input('Which to delete? ([]=rows, other=columns) ','s');
if isempty(flag),
    n = input('How may rows to delete? ([]=1) ');
    if isempty(n),
        n=1;
    else
        assert(n<=Nx-2, 'Out of dimension!');
    end;
    idx = reshape(1:N_kk,Nx,Ny);
    idx = reshape(idx(1:end-n,:),1,[]);
    Nx = Nx-n;
    n_sq_mat(1,kth) = Nx-1;
else
    n = input('How may columns to delete? ([]=1) ');
    if isempty(n),
        n=1;
    else
        assert(n<=Ny-2, 'Out of dimension!');
    end;
    idx = reshape(1:N_kk,Nx,Ny);
    idx = reshape(idx(:,n+1:end),1,[]);
    Ny = Ny-n;
    n_sq_mat(2,kth) = Ny-1;
end;
x_kk = x_kk(:,idx);
y_kk = y_kk(:,idx);
X_kk = X_kk(:,idx);
x_cell{kth} = x_kk;
y_cell{kth} = y_kk;
X_cell{kth} = X_kk;

% plot new axes and corner points on image
figure(6),image(I);
colormap(map);
hold on;
ind_X = Nx*Ny;
ind_XY = Nx;

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
plot([xo;xX]+1,[yo;yX]+1,'g-','linewidth',1);
plot([xo;xY]+1,[yo;yY]+1,'g-','linewidth',1);
plot(y_kk(1,:)+1,y_kk(2,:)+1,'go');
plot(x_kk(1,:)+1,x_kk(2,:)+1,'r+','markersize',10);
text(xXc + delta * uX(1) +1, yXc + delta * uX(2)+1,'X','color','g','Fontsize',14,'HorizontalAlignment','center');
text(xYc + delta * uY(1)+1, yYc + delta * uY(2)+1,'Y','color','g','Fontsize',14,'HorizontalAlignment','center');
text(xo + delta * uO(1) +1, yo + delta * uO(2)+1,'O','color','g','Fontsize',14,'HorizontalAlignment','center');
title(['Check new corners of (camera ' num2str(pp) ', image ' num2str(kk) ')']);
hold off;
disp('Done.');
