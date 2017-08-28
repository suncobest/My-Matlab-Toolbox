% The space size (Unit:mm): (parallel to the initial camera system)
Xmin = 0;
Xmax = 10000;
Ymin = 0;
Ymax = 4000;
Zmin = 0 ;
Zmax = 10000;
center = [Xmin+Xmax; Ymin+Ymax; Zmin+Zmax]/2;

% the 1D calib rig:
% number of points on the calibration rig
np1D = 3;
% np1D = 1+randi(4);

% the 1D coordinates on the rig
rodlen = cumsum(0:np1D-1)*150;
% rodlen = cumsum(0:np1D-1)*100;

flag = input('Add noise to projection or not? ([]=no, other=yes) ','s');
flag = ~isempty(flag);
if flag,
    sigstd = 0.1;   % standard deviation of pixel projection
end;

% motion and interpolation of the stick
gapx = (Xmax-Xmin)/6;
gapy = (Ymax-Ymin)/6;
gapz = (Zmax-Zmin)/6;
n0 = 5;  % number of moving stick
n1 = 100;
n2 = n0*n1;
x0 = reshape([randi(round([Xmin+gapx, Xmax-gapx]),1,n2);
              randi(round([Ymin+gapy, Ymax-gapy]),1,n2);
              randi(round([Zmin+gapz, Zmax-gapz]),1,n2)],n0*3,n1);
% sort x0
for pp=1:n0,
    ii = (pp-1)*3;
    x = x0(ii+1:ii+3,:);
    id = 1;
    index = ones(1,n1);
    restid = true(1,n1);
    restid(id) = 0;
    for kk=2:n1,
        ind = find(restid);
        [~,id] = min(sum(abs(x(:,id)*ones(1,length(ind))-x(:,ind)),1));
        id = ind(id);
        restid(id) = 0;
        index(kk) = id;
    end;
    x0(ii+1:ii+3,:) = x(:,index);
end;

% the origin of rods
s = 20;
n = n1*s;
n_ima = n2*s;
Xori = reshape(permute(reshape(spline_interp3(1:n1,x0,linspace(1,n1,n)),3,n0,n),[1,3,2]),3,n_ima);
q = squad(1:n2,trans_quat_axis(randn(3,n2)),linspace(1,n2,n_ima));

% the direction of rods
Xdir = zeros(3,n_ima);
for kk=1:n_ima,
    R = trans_quat_mat(q(:,kk));
    Xdir(:,kk) = R(:,3);
end;

% all 3D points on rods
Xrod = reshape(permute(repmat(Xori,[1,1,np1D])+Xdir.*reshape(rodlen,[1,1,np1D]), [1,3,2]), 3,[]);

%% cameras:
n_cam = 20;
% n_cam = 5+randi(20);

% intrinsic parameters:
%        fc_mat: focal length of cameras
%        cc_mat: Principal point coordinates
%        alpha_vec: Skew coefficient
%        kc_mat: Distortion coefficients

% imsize = 800+randi(1000,2,n_cam);    % size of CMOS
% fov_angle = 75+randn(1,n_cam)*5;
% fc_mat = imsize(1,:)./tan(pi*fov_angle/360)/2.*[1;1]+randn(2,n_cam)*50;
% cc_mat = (imsize-1)/2+randn(2,n_cam)*50;
% alpha_vec = randi([-1,1],1,n_cam).*rand(1,n_cam)/10;
% kc_mat = [randn(1,n_cam)/10; randn(1,n_cam)/50; randn(2,n_cam)/100; randn(1,n_cam)/1000];

imsize = (randi(1000,2,1)+800)*ones(1,n_cam);
fov_angle = 75+randn*5;
if fov_angle < 30,
    fov_angle = 40;
end;
fc_mat = imsize(1,:)./tan(pi*fov_angle/360)/2.*[1;1]+randn(2,1)*50;
cc_mat = (imsize-1)/2+randn(2,1)*50;
alpha_vec = randi([-1,1]).*rand/10*ones(1,n_cam);
% alpha_vec = zeros(1,n_cam);
kc_mat = zeros(5, n_cam);

% extrinsic parameters:
% the orientation of all camera system

% theta = -randi([10,30],1,n_cam)*pi/180;
theta = -20*pi/180*ones(1,n_cam);
delta = 2*pi/n_cam;
psi = 0:delta:2*pi-delta;
phi = 2*pi*rand(1,n_cam);
Omcw = NaN(3,n_cam);
for pp=1:n_cam,
    Omcw(:,pp) = -rodrigues(trans_euler_mat([psi(pp);theta(pp);phi(pp)],'YXZ'));
end;
% Omcw = randn(3,n_cam);

% handedness of all cameras relative to world coordinate system
hand_list = ones(1,n_cam);
% hand_list = sign(randn(1,n_cam));

% aims of all cameras
x0 = [center(1); center(2)+(Ymax-Ymin)/4; center(3)]*ones(1,n_cam);

% x0 = [center(1)+randi(round([-gapx, gapx]/2),1,n_cam);
%       center(2)+randi(round([-gapy, gapy]/2),1,n_cam);
%       center(3)+randi(round([-gapz, gapz]/2),1,n_cam)];

% the boundary of environment: A*(x;y;z;1)=0, where [x;y;z] = x0-l*diag([1,1,h])*Rt(:,3);
A = [1, 0, 0, -Xmin;
     1, 0, 0, -Xmax;
     0, 1, 0, -Ymin;
     0, 1, 0, -Ymax;
     0, 0, 1, -Zmin;
     0, 0, 1, -Zmax];

% the origin of reference frame in all camera system
Tcw = zeros(3,n_cam);
for pp=1:n_cam,
    Rt = rodrigues(Omcw(:,pp))';
    if hand_list(pp)~=1,
        Rt(3,:) = -Rt(3,:);
    end;
    l = (A(:,1:3)*x0(:,pp)+A(:,4))./(A(:,1:3)*Rt(:,3));
    % the first intersection with boundary
    Tcw(:,pp) = Rt'*(min(l(l>0))*Rt(:,3)-x0(:,pp));
end;


%% 2D projection: generation of 2D points
n_view = n_ima * n_cam;
Np = np1D*n_ima;
active_imgviews = false(n_cam,n_ima);
x_cell = cell(1,n_view);
for pp = 1:n_cam,
    nx = imsize(1,pp);
    ny = imsize(2,pp);
    fc = fc_mat(:,pp);
    cc = cc_mat(:,pp);
    kc = kc_mat(:,pp);
    alpha_c = alpha_vec(pp);
    handkk = hand_list(pp);
    omwkk = Omcw(:,pp);
    Twkk = Tcw(:,pp);
    xx = project_points_mirror2(Xrod,omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
    if flag,
        xx = xx+randn(2,Np)*sigstd;       % add noise
    end;
    mask = reshape(xx(1,:)>-.5 & xx(1,:)<nx-.5 & xx(2,:)>-.5 & xx(2,:)<ny-.5, np1D,n_ima);
    id = all(mask,1);
    active_imgviews(pp,:) = id;
    for kk = find(id),
        kth = (kk-1)*n_cam+pp;
        jj = (kk-1)*np1D;
        x_cell{kth} = xx(:,jj+1:jj+np1D);
    end;
end;

% set bad views of calibration stick inactive (if the pixel distance is less than 20)
delta = 20;
ind_active_views = find(active_imgviews(:)');
for kth = ind_active_views,
    x_kk = x_cell{kth};
    if min(sum(abs(diff(x_kk,[],2)),1))<delta,
        x_cell{kth} = [];
        active_imgviews(kth) = 0;
    end;
end;

active_images = sum(active_imgviews,1)>1;
active_imgviews(:,~active_images) = 0;    % deactivate images with only one camera available
ind_active = find(active_images);
fprintf(1,'Saving generated parameters for one dimensional calibraion!\n');
string_save = ['save multicam_simu_data active_imgviews active_images ind_active Np n_ima ' ...
                   'n_cam n_view np1D rodlen Xmin Xmax Ymin Ymax Zmin Zmax Xrod Xori Xdir '...
                   'imsize fc_mat cc_mat kc_mat alpha_vec hand_list Omcw Tcw x_cell'];
eval(string_save);
write_simu_data;
disp('done.');
