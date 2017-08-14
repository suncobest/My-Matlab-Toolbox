% the 1D calib rig:
% number of points on the calibration rig
np1D = 3;

% the 1D coordinates on the rig
lamda = cumsum(0:np1D-1)*150;

% layout of the blocks
% block_oij = [0, 0; 1, 0; 0, 1];   % origin of blocks (unit: block)
block_oij = [0, 0];

n_block = size(block_oij,1);   % number of blocks

nps = 2; % number of cameras on each side of a block
ncpb = nps*4;
n_cam = n_block*ncpb;

flag = input('Add noise to projection or not? ([]=no, other=yes) ','s');
flag = ~isempty(flag);
if flag,
    sigstd = 0.1;   % standard deviation of pixel projection
end;

% imsize = 800+randi(1000,2,n_cam);    % size of CMOS
% fov_angle = 60+randn(1,n_cam)*10;
% fc_mat = imsize(1,:)./tan(pi*fov_angle/360)/2.*[1;1]+randn(2,n_cam)*50;
% cc_mat = (imsize-1)/2+randn(2,n_cam)*50;
% alpha_vec = (randi(2,1,n_cam)-1).*rand(1,n_cam)/10;
% kc_mat = [randn(1,n_cam)/10; randn(1,n_cam)/50; randn(2,n_cam)/100; randn(1,n_cam)/1000];

imsize = (randi(1000,2,1)+800)*ones(1,n_cam);
fov_angle = 60+randn*10;
if fov_angle < 30,
    fov_angle = 40;
end;
fc_mat = imsize(1,:)./tan(pi*fov_angle/360)/2.*[1;1]+randn(2,1)*50;
cc_mat = (imsize-1)/2+randn(2,1)*50;
alpha_vec = (randi(2,1)-1).*rand/10*ones(1,n_cam);
% alpha_vec = zeros(1,n_cam);
kc_mat = zeros(5, n_cam);

nfpb = 2000;  % number of frames in each block
n_ima = nfpb*n_block;
Xdir = randn(3,n_ima);   % the direction of rods
Xdir = Xdir./(ones(3,1)*sqrt(sum(Xdir.^2,1)));
Xori = NaN(3,n_ima);   % the origin of rods

% The block space size (Unit:mm): (parallel to the initial camera system)
Xlen = 7000;
Ylen = 3000;
Zlen = 7000;
Omcw = NaN(3,n_cam);
Tcw = Omcw;
% handedness of all cameras relative to world coordinate system
hand_list = ones(1,n_cam);
% hand_list = sign(randn(1,n_cam));

for count = 1:n_block,
    % origin of current block
    x0 = [Xlen*block_oij(count,1); 0; Zlen*block_oij(count,2)];
    Xori(:,(count-1)*nfpb+1:count*nfpb) = [Xlen*rand(1,nfpb); Ylen*rand(1,nfpb); Zlen*rand(1,nfpb)] + x0;

    % the origin of cameras in reference system
    Tc = [Xlen*(1:nps)/(nps+1)+randn(1,nps)*100, Xlen*ones(1,nps), Xlen*(nps:-1:1)/(nps+1)+randn(1,nps)*100, zeros(1,nps);
          zeros(1,ncpb);
          zeros(1,nps), Zlen*(1:nps)/(nps+1)+randn(1,nps)*100, Zlen*ones(1,nps), Zlen*(nps:-1:1)/(nps+1)+randn(1,nps)*100] + x0;

    % the orientation of all camera system
    yaw = [zeros(1,nps), -pi/2*ones(1,nps), pi*ones(1,nps), pi/2*ones(1,nps)]+(-0+0*rand(1,ncpb))*pi/180;
    pitch = -(25+0*rand(1,ncpb))*pi/180;
    roll = (-15+30*rand(1,ncpb))*pi/180;
    % roll = randn(1,ncpb)/10;    % std = 5.7 degree
    ii = (count-1)*ncpb;
    for pp=1:ncpb,
        jj = ii+pp;
        R = trans_euler_mat([yaw(pp);pitch(pp);roll(pp)],'YXZ')';
        Omcw(:,jj) = rodrigues(R);
        if hand_list(jj)~=1,
            R(:,3) = -R(:,3);
        end;
        % the origin of reference frame in all cameras
        Tcw(:,jj) = -R*Tc(:,pp);
    end;
end;

% all 3D points on rods
Xrod = reshape(permute(Xori+Xdir.*reshape(lamda,[1,1,np1D]), [1,3,2]), 3,[]);

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
    mask = reshape(xx(1,:)>-.5 & xx(1,:)<nx-.5 & xx(2,:)>-.5 & xx(2,:)<ny-.5, np1D, n_ima);
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
                   'n_cam n_view np1D lamda Xlen Ylen Zlen Xrod Xori Xdir hand_list Omcw Tcw '...
                   'imsize fc_mat cc_mat kc_mat alpha_vec x_cell'];
eval(string_save);
write_simu_data;
disp('done.');
