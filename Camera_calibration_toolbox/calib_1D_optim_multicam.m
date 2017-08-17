% calib_1D_optim_multicam
%
% see calib_1D_optim_mirror.
%
% Main calibration function. Computes the intrinsic and extrinsic parameters.
% Runs as a script.
%
% INPUT: x_cell: Feature locations on the images;
%        rodlen: the 1D coordinates of points on the rod;
%
% OUTPUT: fc_mat: focal length of cameras;
%        cc_mat: Principal point coordinates;
%        alpha_vec: Skew coefficient;
%        kc_mat: Distortion coefficients;
%        Omcc: axis angle rotation vectors of camera 1 to other cameras;
%        Tcc: 3D translation vectors from camera 1 to other cameras;
%        Xorg: the origin of reconstructed sticks in camera 1;
%        Xdir: the direction of reconstructed sticks;
%        Xrod: the 3D coordinates of points on reconstructed sticks;
%
% Method: Minimizes the pixel reprojection error in the least squares sense over the intrinsic
%        camera parameters, and the extrinsic parameters (3D locations of the grids in space)
%
% VERY IMPORTANT: This function works for 1D calibration rig.

%  By ZPF @ZVR, 2017-7-27

if ~exist('x_cell','var') || ~exist('imsize','var'),
    if exist('multicam_simu_data.mat','file')==2,
        fprintf(1,'\nSimulation data detected! Do you want to load it?\n');
        flag = input('Load the simulation data or not? ([]=yes, other=no) ','s');
        if isempty(flag),
            load('multicam_simu_data.mat');
            % clear Xrod Xori Xdir fc_mat cc_mat kc_mat alpha_vec Omcw Tcw;
        end;
    else
        fprintf(1,'\nThere is no data required for calibration!\n');
        return;
    end;
end;

for kk = ind_active,
    if all(active_imgviews(:,kk)==0),
        fprintf(1,'No camera captured points in frame %d, The image is now set inactive!\n',kk);
    end;
end;
active_images = sum(active_imgviews,1)>1;
ind_active = find(active_images);
active_imgviews(:,~active_images) = 0;

if ~exist('MaxIter','var'),
    MaxIter = 30; % Maximum number of iterations in the main LM algorithm
end;

if ~exist('check_cond','var'),
    check_cond = 1;
end;

if ~exist('desactivated_images','var'),
    desactivated_images = [];
end;

if ~exist('est_alpha_vec','var'),
    flag = 1;
    fprintf(1,'\nVector of pixel skew estimation indicator do not exist!\n');
else
    est_alpha_vec = est_alpha_vec(:)';
    flag = length(est_alpha_vec)~=n_cam;
    if flag,
        fprintf(1,'\nDimension of est_alpha_vec do not match number of cameras! Please input again!\n');
    end;
end;
while flag,
    fprintf(1,'\nDo you want to estimate the pixel skew of all %d cameras?\nSet to zero if you don''t!\n',n_cam);
    % by default do not estimate skew
    est_alpha_vec = input(['est_alpha_vec = ([] = [' num2str(zeros(1,n_cam)) '])']);
    if isempty(est_alpha_vec),
        est_alpha_vec = zeros(1,n_cam);
        flag = 0;
    else
        est_alpha_vec = est_alpha_vec(:)';
        flag = length(est_alpha_vec)~=n_cam;
        if flag,
            fprintf(1,'\nDimension of est_alpha_vec do not match number of cameras! Please input again!\n');
        end;
    end;
end;

if ~exist('est_aspect_ratio_vec','var'),
    flag = 1;
    fprintf(1,'\nVector of aspect ratio estimation indicator do not exist!\n');
else
    est_aspect_ratio_vec = est_aspect_ratio_vec(:)';
    flag = length(est_aspect_ratio_vec)~=n_cam;
    if flag,
        fprintf(1,'\nDimension of est_aspect_ratio_vec do not match number of cameras! Please input again!\n');
    end;
end;
while flag,
    fprintf(1,'\nDo you want to estimate the aspect ratio of all %d cameras?\nSet to zero if you don''t!\n',n_cam);
    % by default estimate aspect ratio
    est_aspect_ratio_vec = input(['est_aspect_ratio_vec = ([] = [' num2str(ones(1,n_cam)) '])']);
    if isempty(est_aspect_ratio_vec),
        est_aspect_ratio_vec = ones(1,n_cam);
        flag = 0;
    else
        est_aspect_ratio_vec = est_aspect_ratio_vec(:)';
        flag = length(est_aspect_ratio_vec)~=n_cam;
        if flag,
            fprintf(1,'\nDimension of est_aspect_ratio_vec do not match number of cameras! Please input again!\n');
        end;
    end;
end;

if ~exist('center_optim_vec','var'),
    flag = 1;
    fprintf(1,'\nVector of principal point estimation indicator do not exist!\n');
else
    center_optim_vec = center_optim_vec(:)';
    flag = length(center_optim_vec)~=n_cam;
    if flag,
        fprintf(1,'\nDimension of center_optim_vec do not match number of cameras! Please input again!\n');
    end;
end;
while flag,
    fprintf(1,'\nDo you want to estimate the principal point of all %d cameras?\nSet to zero if you don''t!\n',n_cam);
    % by default estimate principal points
    center_optim_vec = input(['center_optim_vec = ([] = [' num2str(ones(1,n_cam)) '])']);
    if isempty(center_optim_vec),
        center_optim_vec = ones(1,n_cam);
        flag = 0;
    else
        center_optim_vec = center_optim_vec(:)';
        flag = length(center_optim_vec)~=n_cam;
        if flag,
            fprintf(1,'\nDimension of center_optim_vec do not match number of cameras! Please input again!\n');
        end;
    end;
end;

if ~exist('est_fc_mat','var'),
    flag = 1;
    fprintf(1,'\nMatrix of focal length estimation indicator do not exist!\n');
else
    flag = ~isequal(size(est_fc_mat),[2,n_cam]);
    if flag,
        fprintf(1,'\nDimension of est_fc_mat do not match number of cameras! Please input again!\n');
    end;
end;
if flag,
    est_fc_mat = ones(2,n_cam);
    flag = input('Estimate all cameras'' focal length or not? ([]=yes, other=no) ','s');
    if ~isempty(flag),
        for pp = 1:n_cam,
            flag = 1;
            while flag,
                fprintf(1,'\nDo you want to estimate the focal length of camera %d?\nSet to zero if you don''t!\n',pp);
                % by default estimate focal length
                est_fc = input(['est_fc [fc1; fc2] of camera ' num2str(pp) ': ([] = [1;1])])']);
                if isempty(est_fc),
                    est_fc = [1;1];
                    flag = 0;
                else
                    est_fc = est_fc(:);
                    flag = length(est_fc)~=2;
                    if flag,
                        fprintf(1,'\nUnexpected dimension of est_fc! Please input again!\n');
                    end;
                end;
            end;
            est_fc_mat(:,pp) = est_fc;
        end;
    end;
end;

if ~exist('est_dist_mat','var'),
    flag = 1;
    fprintf(1,'\nMatrix of distortion coefficients estimation indicator do not exist!\n');
else
    flag = ~isequal(size(est_dist_mat),[5,n_cam]);
    if flag,
        fprintf(1,'\nDimension of est_dist_mat do not match number of cameras! Please input again!\n');
    end;
end;
if flag,
    est_dist_mat = zeros(5,n_cam);
    flag = input('Estimate lens distortion or not? ([]=no, other=yes) ','s');
    if ~isempty(flag),
        for pp = 1:n_cam,
            fprintf(1,'\nDo you want to estimate the distortion coefficients of camera %d?\nSet to zero if you don''t!\n',pp);
            est_dist = input(['est_dist [k1; k2; k3; k4; k5] of camera ' num2str(pp) ': ([] = [1;1;0;0;0])])']);
            if isempty(est_dist),
                est_dist = [1;1;0;0;0];
            else
                est_dist = est_dist(:);
                n = length(est_dist);
                if n<5,
                    est_dist = [est_dist;zeros(5-n,1)];
                elseif n>5,
                    est_dist = est_dist(1:5);
                end;
            end;
            fprintf(1,'\nest_dist of camera %d: [%d; %d; %d; %d; %d].\n',pp, est_dist);
            est_dist_mat(:,pp) = est_dist;
        end;
    end;
end;

% Little fix in case of stupid values in the binary variables:
center_optim_vec = double(~~center_optim_vec);
est_alpha_vec = double(~~est_alpha_vec);
est_aspect_ratio_vec = double(~~est_aspect_ratio_vec);
est_dist_mat = double(~~est_dist_mat);
est_fc_mat = double(~~est_fc_mat);

fprintf(1,'\n');

if ~exist('fc_mat','var'),
    flag = input('The calibration rod was only under rotation or not? ([]=no, other=yes) ','s');
    if isempty(flag),
        FOV_angle = 90; %field of view in degrees: for 135 camera, 90 degree of FOV is about 18 mm focal length。
        fprintf(1,'Initialization of the focal length with FOV of %3.1f degrees.\n\n',FOV_angle);
        fc_mat = ones(2,1)*(imsize(1,:)/2)/tan(pi*FOV_angle/360);    % FOV_angle=2*atan(nx/(2*fc))
    else
        fprintf(1,'\nInitialization of the intrinsic parameters using Zhang Zhengyou''s algorithm.\n');
        intrinsic_1D_rot_multicam;
    end;
end;

% Initialization of the intrinsic parameters
if ~exist('cc_mat','var'),
    fprintf(1,'\nInitialization of the principal point at the center of each camera view.\n');
    cc_mat = (imsize-1)/2;     % cc = [(nx-1)/2;(ny-1)/2];
end;

if ~exist('kc_mat','var'),
    fprintf(1,'\nInitialization of all camera image distortion to zero.\n');
    kc_mat = zeros(5,n_cam);
else
    [m,n] = size(kc_mat);
    if n~=n_cam,
        fprintf(1,'\nDimension of kc_mat do not match number of cameras! Re-initialization kc_mat to zero!\n');
        kc_mat = zeros(5,n_cam);
    end;
    if m<5,
        fprintf(1,'\nRadial distortion coefficient up to the 6th degree.\n');
        kc_mat = [kc_mat;zeros(5-m,n_cam)];
    elseif m>5,
        kc_mat = kc_mat(1:5,:);
    end;
end;

if ~exist('alpha_vec','var'),
    fprintf(1,'Initialization of all camera image skew to zero.\n');
    alpha_vec = zeros(1,n_cam);
end;

% no principal point estimation, no skew estimation
est_alpha_vec(center_optim_vec==0) = 0;

% alpha = 0 when skew is not estimated
alpha_vec(est_alpha_vec==0) = 0;

% set to zero the distortion coefficients that are not estimated
kc_mat = kc_mat .* est_dist_mat;

% Conditioning threshold for view rejection
thresh_cond = 1e4;
% threshold to terminate the main LM iteration
gradeps = eps*1e8;


%%% calculate the shortest path from all cameras to selected camera.
A_cam = zeros(n_cam);   % adjacency matrix of all camera nodes
for pp = 1:n_cam,
    for kk = pp+1:n_cam,
        tmp = 1/sum(all(active_imgviews([pp,kk],:),1));   % weight of path
        A_cam(pp,kk) = tmp;
        A_cam(kk,pp) = tmp;
    end;
end;
[costs, paths] = dijkstra(A_cam);
[~, idm] = min(sum(costs, 2));
pathm = paths(idm,:);

%%% Initialization of the camera parameters for global minimization:
%%% Computes pairwise extrinsic parameters for all cameras
nc = 10;  % threshold number of common images
n_view = n_ima * n_cam;
Omcc = NaN(3, n_cam);  % transformation from the main camera to all cameras
Tcc = NaN(3, n_cam);
Omcc(:,idm) = 0;
Tcc(:,idm) = 0;
active_cam = true(1,n_cam);
handcc = hand_list(idm)*hand_list;

for pp = [1:idm-1, idm+1:n_cam],
    if active_cam(pp),
        camlist = pathm{pp};
        for ii = 1:length(camlist)-1,
            id = camlist([ii,ii+1]);  % index of camera pairs
            if ~isnan(Omcc(1,id(2))),
                % relation of the two cameras is set
                continue;
            end;
            ind = all(active_imgviews(id,:),1);  % logical index of common images
            flag = sum(ind,2) >= nc;
            if flag,
                % compute pairwise camera relation: ||T2||=1
                fc = fc_mat(:,id);
                cc = cc_mat(:,id);
                kc = kc_mat(:,id);
                alpha_c = alpha_vec(id);
                handkk = handcc(id(1))*handcc(id(2));
                kk = find(ind);
                xx = cell2mat(x_cell(repmat((kk-1)*n_cam,[1,1,2])+reshape(id,[1,1,2])));
                [om2,T2] = compute_Rt_pair(xx,fc,cc,kc,alpha_c,handkk);

                % triangulation and  determine the scale factor
                XX = compute_structure2(xx,[zeros(3,1),om2],[zeros(3,1),T2],[1,handkk],fc,cc,kc,alpha_c);
                XX = reshape(XX,[3,np1D,length(kk)]);
                ind = all(all(~isnan(XX),1),2);
                Xlen = permute(sqrt(sum(diff(XX(:,:,ind),[],2).^2,1)),[3,2,1]);  % length of rod (with ||T2||=1)
                s = mean(diff(rodlen)./mean(Xlen,1));

                % chain the Euclidean motion
                [om,T] = compose_motion2(Omcc(:,id(1)),Tcc(:,id(1)),om2,T2*s,handkk);
                Omcc(:,id(2)) = om;
                Tcc(:,id(2)) = T;
            else
                active_cam(camlist(ii+1:end)) = 0;
                break;
            end;
        end;
    end;
end;

active_imgviews(~active_cam,:) = 0;
active_images = sum(active_imgviews,1)>1;
ind_active = find(active_images);
active_imgviews(:,~active_images) = 0;
active_cam = any(active_imgviews,2)';
ind_cam = find(active_cam);
nc = length(ind_cam);

% reconstruct the calibration rod in the main camera:
npts = n_ima*np1D;
ind_active_views = find(active_imgviews(:)');
xx = NaN(2,np1D,n_cam,n_ima);
xx(:,:,ind_active_views) = reshape(cell2mat(x_cell(ind_active_views)),2,np1D,[]);
xx = reshape(permute(xx,[1,2,4,3]),[2,npts,n_cam]);
XX = compute_structure2(xx,Omcc,Tcc,handcc,fc_mat,cc_mat,kc_mat,alpha_vec);
Xo = XX(:,1:np1D:end);
Xn = XX(:,np1D:np1D:end)-Xo;
Xn = Xn./(ones(3,1)*sqrt(sum(Xn.^2,1)));
thph = cartesian2spherical(Xn);
thph = thph(2:3,:);

if exist('Omcw','var'),
    Om2 = NaN(3,n_cam);
    T2 = Om2;
    err = zeros(6,nc);
    for kk = 1:nc,
        pp = ind_cam(kk);
        [om,T] = compose_motion2(Omcw(:,idm),Tcw(:,idm),Omcc(:,pp),Tcc(:,pp),handcc(pp));
        Om2(:,pp) = om;
        T2(:,pp) = T;
        err(1:3,kk) = om-Omcw(:,pp);
        err(4:6,kk) = (T-Tcw(:,pp))/norm(Tcw(:,pp));
    end;
    if exist('Xrod','var'),
        ind = reshape(repmat(active_images,[np1D,1]),1,npts);
        errX = XX(:,ind)-rigid_refmotion(Xrod(:,ind),Omcw(:,idm),Tcw(:,idm),hand_list(idm));
    end;
end;


%% ------------------------------------------ Main Optimization:

fprintf(1,'\nMain calibration optimization procedure - Bundle adjustment with %d cameras\n', n_cam);
fprintf(1,'Sparse Levenberg-Marquardt iterations: ');

%%% Initialization of the global parameter vector:
ncam16 = 16*n_cam;
nima5 = 5*n_ima;
nima3 = 3*n_ima;
nima2 = 2*n_ima;
npts3 = 3*npts;
intr_update = reshape([fc_mat; cc_mat; alpha_vec; kc_mat; Omcc; Tcc],ncam16,1);
extr_update = reshape([Xo; thph],nima5,1);
init_param = [intr_update; extr_update];

% initial error before bundle adjustment
intr_param = init_param(1:ncam16);
extr_param = init_param(ncam16+1:end);
Xp = reshape(extr_param,5,n_ima);
Xp = gen_1D_points(Xp(1:3,:),Xp(4:5,:),rodlen);
ex = []; % Global error vector
for pp = ind_cam,
    ii = (pp-1)*16;
    % load camera parameters
    fc = intr_param(ii+1 : ii+2);
    cc = intr_param(ii+3 : ii+4);
    alpha_c = intr_param(ii+5);
    kc = intr_param(ii+6 : ii+10);
    omwkk = intr_param(ii+11 : ii+13);
    Twkk = intr_param(ii+14 : ii+16);
    handkk = handcc(pp);
    active_view = active_imgviews(pp,:);
    kth = (find(active_view)-1)*n_cam+pp;
    x_kk = cell2mat(x_cell(kth));
    ind = reshape(repmat(active_view,[np1D,1]),1,npts);
    x = project_points_mirror2(Xp(:,ind),omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
    ex_kk = x_kk - x;
    ex = [ex, ex_kk];
end;
err_std0 = std(ex,0,2);

% return;

% The following vector helps to select the variables to update:
selected_invar = [est_fc_mat; ones(2,1)*center_optim_vec; est_alpha_vec; est_dist_mat; active_cam(ones(6,1),:)];
selected_invar(2,:) = selected_invar(2,:).*(est_aspect_ratio_vec | ~est_fc_mat(1,:));
selected_invar(11:16,idm) = 0;
selected_exvar = reshape(active_images(ones(5,1),:),1,nima5);
ind_va = find(selected_invar);
ind_vb = find(selected_exvar);

lamda = 0.001; % set an initial value of the damping factor for the LM
updateJ = 1;
ex = ex(:);
ex2_init = dot(ex,ex);
ex2 = ex2_init;
for iter = 1:MaxIter,
    fprintf(1,'%d...',iter);
    if mod(iter,20)==0,
        fprintf(1,'\n');
    end;
    if updateJ,
        % JJ2 = JJ'*JJ = [U, W; W', V]
        U = sparse([],[],[],ncam16,ncam16,16*ncam16);
        V = sparse([],[],[],nima5,nima5,5*nima5);
        W = sparse([],[],[],ncam16,nima5,ncam16*nima5);
        ea = zeros(ncam16,1);        % A'*ex
        eb = zeros(nima5,1);         % B'*ex

        % generate 1D points on rod
        Xp = reshape(extr_param,5,n_ima);
        Xo = Xp(1:3,:);
        thph = Xp(4:5,:);
        [Xp,dXpdXo,dXpdtp] = gen_1D_points(Xo,thph,rodlen);
        for pp = ind_cam,
            ii = (pp-1)*16;
            % load camera parameters
            fc = intr_param(ii+1 : ii+2);
            cc = intr_param(ii+3 : ii+4);
            alpha_c = intr_param(ii+5);
            kc = intr_param(ii+6 : ii+10);
            omwkk = intr_param(ii+11 : ii+13);
            Twkk = intr_param(ii+14 : ii+16);
            handkk = handcc(pp);

            % load pixel points
            active_view = active_imgviews(pp,:);
            kth = (find(active_view)-1)*n_cam+pp;
            x_kk = cell2mat(x_cell(kth));

            ind = reshape(repmat(active_view,[np1D,1]),1,npts);
            if est_aspect_ratio_vec(pp),
                [x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha,dxdXp] = project_points_mirror2(Xp(:,ind),omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
            else
                [x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha,dxdXp] = project_points_mirror2(Xp(:,ind),omwkk,Twkk,handkk,fc(1),cc,kc,alpha_c);
                dxdf = repmat(dxdf,[1 2]);
            end;
            ex_kk = x_kk - x;
            Akk = [dxdf dxdc dxdalpha dxdk dxdom dxdT];

            na = sum(active_view);
            nb = 2*na*np1D;
            idx = reshape(repmat(ind,[3,1]),1,npts3);
            id3 = reshape(repmat(active_view,[3,1]),1,nima3);
            id2 = reshape(repmat(active_view,[2,1]),1,nima2);
            id5 = reshape(repmat(active_view,[5,1]),1,nima5);
            Bkk = reshape([reshape(dxdXp*dXpdXo(idx,id3),[nb,3,na]), reshape(dxdXp*dXpdtp(idx,id2),[nb,2,na])], nb,5*na);

            U(ii+1 : ii+16, ii+1 : ii+16) = Akk'*Akk;
            ea(ii+1 : ii+16) = Akk'*ex_kk(:);
            V(id5, id5) = V(id5, id5) + Bkk'*Bkk;
            W(ii+1 : ii+16, id5) = Akk'*Bkk;
            eb(id5) = eb(id5) + Bkk'*ex_kk(:);
        end;
        U = U(ind_va,ind_va);
        V = V(ind_vb,ind_vb);
        W = W(ind_va,ind_vb);
        ea = ea(ind_va);
        eb = eb(ind_vb);
    end;

    U_lm = U + diag(lamda*diag(U));  % U + lamda*speye(size(U));
    V_lm = V + diag(lamda*diag(V));  %  V + lamda*speye(size(V));
    Y = W/V_lm;

    intr_innov = (U_lm-Y*W')\(ea-Y*eb);              % da
    extr_innov = V_lm\(eb-W'*intr_innov);            % db

    intr_update(ind_va) = intr_param(ind_va) + intr_innov;     % updated parameters
    extr_update(ind_vb) = extr_param(ind_vb) + extr_innov;

    % generate 1D points on rod
    Xp = reshape(extr_update,5,n_ima);
    Xp = gen_1D_points(Xp(1:3,:),Xp(4:5,:),rodlen);
    % compute reprojection error vector
    ex = [];
    for pp = ind_cam,
        ii = (pp-1)*16;
        % New intrinsic parameters:
        if ~est_aspect_ratio_vec(pp) && all(est_fc_mat(:,pp)),
            intr_update(ii+2) = intr_update(ii+1);
        end;
        fc = intr_update(ii+1 : ii+2);
        cc = intr_update(ii+3 : ii+4);
        alpha_c = intr_update(ii+5);
        kc = intr_update(ii+6 : ii+10);
        handkk = handcc(pp);
        omwkk = intr_update(ii+11 : ii+13);
        Twkk = intr_update(ii+14 : ii+16);
        if center_optim_vec(pp) && (cc(1)<-.5 || cc(1)>imsize(1,pp)-.5 || cc(2)<-.5 || cc(2)>imsize(2,pp)-.5),
            center_optim_vec(pp) = 0;
            fprintf(1,['Warning: it appears that the principal point of camera %d'...
                           'cannot be estimated.\nPlease run this code again!\n'], pp);
            return;
        end;
        active_view = active_imgviews(pp,:);
        kth = (find(active_view)-1)*n_cam+pp;
        x_kk = cell2mat(x_cell(kth));
        ind = reshape(repmat(active_view,[np1D,1]),1,npts);
        x = project_points_mirror2(Xp(:,ind),omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
        ex_kk = x_kk - x;
        ex = [ex, ex_kk];
    end;
    ex = ex(:);
    ex2_update = dot(ex,ex);
    if ex2_update < ex2,
        lamda = lamda/10;
        intr_param = intr_update;
        extr_param = extr_update;
        ex2 = ex2_update;
        updateJ=1;
    else
        lamda = lamda*10;
        updateJ=0;
    end;
    grad = -[ea;eb];
    if norm(grad) < gradeps,
        break;
    end;
end;

% return;

%%%--------------------------- Computation of the error of estimation:

fprintf(1,'\nEstimation of uncertainties...');

% Extraction of the paramters for computing the reprojection error:
intr_param = reshape(intr_param, 16, n_cam);
fc_mat = intr_param(1:2,:);
cc_mat = intr_param(3:4,:);
alpha_vec = intr_param(5,:);
kc_mat = intr_param(6:10,:);
Omcc = extr_param(11:13,:);
Tcc = extr_param(14:16,:);

extr_param = reshape(extr_param, 5, n_ima);
Xo = extr_param(1:3,:);
thph = extr_param(4:5,:);
Xp = gen_1D_points(Xo,thph,rodlen);

% compute the reprojected errors
y_cam = cell(1, n_cam);  % Reprojected points
ex_cam = y_cam;          % Reprojected error
err_cam = zeros(2,n_cam);
ex = [];
for pp = ind_cam,
    fc = fc_mat(:,pp);
    cc = cc_mat(:,pp);
    alpha_c = alpha_vec(pp);
    kc = kc_mat(:,pp);
    handkk = handcc(pp);
    active_view = active_imgviews(pp,:);
    kth = (find(active_view)-1)*n_cam+pp;
    x_kk = cell2mat(x_cell(kth));
    ind = reshape(repmat(active_view,[np1D,1]),1,npts);
    y_kk = project_points_mirror2(Xp(:,ind),omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
    ex_kk = x_kk - y_kk;
    ex = [ex, ex_kk];
    y_cam{pp} = y_kk;
    ex_cam{pp} = ex_kk;
    err_cam(:,pp) = std(ex_kk,0,2);
end;

err_std = std(ex,0,2);
sigma_x = std(ex(:));

% Compute the the standard deviation of parameters from std(ex):
% ex = X-f(P),  cov(param,param) = inv((JJ'* inv(cov(ex,ex))* JJ))
JJ2 = [U, W; W', V];    % JJ2 = JJ'*JJ
JJ2_inv = speye(size(JJ2))/JJ2;

% Extraction of the final intrinsic paramaters:
param_error = zeros(ncam16,1);           % take value as perfect if not optimized (no error).
nva = length(ind_va);
param_error(ind_va) = 3*sqrt(abs(full(diag(JJ2_inv(1:nva,1:nva)))))*sigma_x;     % 3 sigma principle
param_error = reshape(param_error,16,n_cam);
ind = find(~est_aspect_ratio_vec & all(est_fc_mat,1));
param_error(2,ind) = param_error(1,ind);

fc_mat_error = param_error(1:2,:);
cc_mat_error = param_error(3:4,:);
alpha_vec_error = param_error(5,:);
kc_mat_error = param_error(6:10,:);
Omcc_error = param_error(11:13,:);
Tcc_error = param_error(14:16,:);

% Extraction of the final extrinsic paramaters:
param_error = zeros(nima5,1);
nvb = length(ind_vb);
param_error(ind_vb) = 3*sqrt(abs(full(diag(JJ2_inv(nva+1:nva+nvb, nva+1:nva+nvb)))))*sigma_x;     % 3 sigma principle
param_error = reshape(param_error,5,n_ima);
Xo_error = param_error(1:3,:);
thph_error = param_error(4:5,:);

fprintf(1,'Done with the uncertainties!\n');

for pp = ind_cam,
    fc = fc_mat(:,pp);
    cc = cc_mat(:,pp);
    kc = kc_mat(:,pp);
    alpha_c = alpha_vec(pp);
    fc_error = fc_mat_error(:,pp);
    cc_error = cc_mat_error(:,pp);
    alpha_c_error =  alpha_vec_error(pp);
    kc_error = kc_mat_error(:,pp);
    fprintf(1,'\n\nCalibration results of camera %d after optimization (with uncertainties):\n\n',pp);
    fprintf(1,'Focal Length:      fc = [%3.5f, %3.5f] ? [%3.5f, %3.5f]\n',[fc;fc_error]);
    fprintf(1,'Principal point:   cc = [%3.5f, %3.5f] ? [%3.5f, %3.5f]\n',[cc;cc_error]);
    fprintf(1,'Skew:         alpha_c = %3.5f ? %3.5f => Skew angle = %3.5f ? %3.5f degrees\n', ...
        [alpha_c; alpha_c_error; 90-atan(alpha_c)*180/pi; atan(alpha_c_error)*180/pi]);
    fprintf(1,'Distortion:          kc = [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f] ? [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f]\n',[kc;kc_error]);
    fprintf(1,'Pixel error:        err = [%3.5f, %3.5f]\n\n',err_cam(:,pp));
    fprintf(1,'Note: The pixel error are approximately three times the standard deviations (for reference).\n\n\n');

    % Some recommendations to the user to reject some of the difficult unkowns....
    alpha_c_min = alpha_c - alpha_c_error/2;
    alpha_c_max = alpha_c + alpha_c_error/2;
    if (alpha_c_min < 0) && (alpha_c_max > 0),
        fprintf(1,'Note: the skew coefficient of camera %d is found to be equal to zero (within its uncertainty).\n', pp);
        fprintf(1,'Setting est_alpha_vec(%d)=0;\n\n',pp);
        est_alpha_vec(pp)=0;
    end;
    kc_min = kc - kc_error/2;
    kc_max = kc + kc_error/2;
    prob_kc = (kc_min < 0) & (kc_max > 0);
    if ~(prob_kc(3) && prob_kc(4)),
        prob_kc(3:4) = [0;0];
    end;
    est_dist = est_dist_mat(:,pp);
    if any(prob_kc),
        est_dist = est_dist & ~prob_kc;
        fprintf(1,'Note: Some distortion coefficients of camera %d are found equal to zero (within their uncertainties).\n',pp);
        fprintf(1,'Setting est_dist_mat(:,%d)=[%d;%d;%d;%d;%d];\n\n',pp,est_dist);
        est_dist_mat(:,pp) = est_dist;
    end;
end;
fprintf(1,['For accurate and stable error estimates, it is recommended to run Calibration twice.\n' ...
    'You may using projection to recompute the conners to improve the quality of your data!\n\n']);


%% Refine the camera parameters if more than one camera is available.
% First,determine the relative orientation and position of each view to the first one.
% Then use intrinsic parameters, relations between cameras, and the extrinsic parameters
% of the first view as the initial value for the refinement.

refine_multicam = 0;
if n_cam>1,
    fprintf(1,['More than one camera detected, do you want to trigger the multicamera refinement now?\n'...
                   'If you have already run the Calibration twice, it is recommended to trigger refinement!\n']);
    refine_multicam = input('Refine the multicamera calibration or not? ([]=yes, other=no) ','s');
    refine_multicam = isempty(refine_multicam);
end;
if ~refine_multicam,
    return;
end;

% Initialization of the extrinsic parameters for multi-view refinement:
hand1 = hand_list(1);
handcc = hand1*hand_list;       % handcc is the handedness of every camera frame wrt camera 1
Omcw = NaN(3, n_ima);            % Omcw and Tcw are the rotation and translation of camera 1 wrt every world frame
Tcw = NaN(3, n_ima);

Omcc = NaN(3, n_cam);           % Omcc and Tcc are the rotation and translation of every camera wrt camera 1
Omcc(:,1) = zeros(3,1);
Tcc = NaN(3, n_cam);
Tcc(:,1) = zeros(3,1);

Omcw_cell = cell(1,n_ima);
Tcw_cell = Omcw_cell;
Omcc_cell = cell(1, n_cam);      % Omcc_cell, Tcc_cell store transformation of every camera wrt camera 1
Tcc_cell = Omcc_cell;
ind_active_views = find(active_imgviews(1,:));

for kk = ind_active_views,
    jj = (kk-1)*n_cam+1;
    omw1 = Omw_mat(:, jj);              %   Rc1= rodrigues(omw1);
    Tw1= Tw_mat(:, jj);
    Rw1t = rodrigues(-omw1);
    Omcw(:, kk) = omw1;
    Tcw(:, kk) = Tw1;
    active_view = find(active_imgviews(2:end, kk)')+1;
    for pp = active_view,
        kth = jj+pp-1;
        omwkk = Omw_mat(:, kth);
        Twkk = Tw_mat(:, kth);
        handkk = handcc(pp);
        [omck, Tck] = compose_motion2(-omw1,-Rw1t*Tw1,omwkk,Twkk,handkk);  %  inverse composition
        Omcc_cell{pp} = [Omcc_cell{pp}, omck];
        Tcc_cell{pp} = [Tcc_cell{pp}, Tck];
    end;
end;

% relative deviation checking of close rotation and translation of each camera wrt the 1st one.
bigeps = 1e-2;
extrinsic_deviation = zeros(1,n_cam);
for pp=2:n_cam,
    omck = Omcc_cell{pp};
    Tck = Tcc_cell{pp};
    n = size(omck,2);
    switch n,
        case 0,
            fprintf(1,'No transformation of camera %d wrt camera 1 found!\n',pp);
            active_imgviews(pp,:) = 0;
            continue;
        case 1,
            Omcc(:, pp) = omck;
            Tcc(:, pp) = Tck;
        otherwise,
            % transform axis angle to quaternion to compute mean rotation
            Qck = trans_quat_axis(omck);
            Qc1 = Qck(:,1);
            for i = 2:n,
                Qcki = Qck(:,i);
                if sum(Qcki.*Qc1)<0,
                    Qck(:,i) = -Qcki;
                end;
            end;
            Q = quatmean(Qck,ones(1,n)/n);        % mean of quaternions
            if sum(Q.*Qc1)<0,
                Q = -Q;
            end;
            Tcc(:, pp) = mean(Tck,2);        % mean of translation
            Omcc(:, pp) = trans_quat_axis(Q);
            % relative deviation of orientation and position of cameras
            deviation = max(max(sqrt(sum((Qck-Q(:,ones(1,n))).^2)),[],2), max(sqrt(sum((Tck-Tcc(:,pp)*ones(1,n)).^2)),[],2)/norm(Tcc(:,pp)));
            fprintf(1,'For camera %d, relative extrinsic deviation = %1.0e;\n', pp, deviation);
            if deviation>bigeps,
                fprintf(1,'The value is larger than %1.0e, please stop and check variables ''Qck'' and ''Tck''.\n', bigeps);
                flag = input('Abort calibration or not? ([]=yes, other=no) ','s');
                if isempty(flag),
                    fprintf(1,'Program terminated!\n');
                    return;
                end;
            end;
            extrinsic_deviation(pp) = deviation;
    end;
end;
active_images = any(active_imgviews,1);
ind_active = find(active_images);

% If Xk=Rk*Hk*Xw+Tk, Xk=Rk1*Hk1*X1+Tk1, where Hk1=Hk*H1
% then X1=Hk1*Rk1'*Rk*Hk1*H1*Xw+Hk1*Rk1'*(Tk-Tk1)
% Given om=[n1;n2;n3], Hzkk = diag([1,1, -1]), if Rt=Hzkk*rodrigues(om)*Hzkk,
% then we have omt = rodrigues(Rt) = [-n1; -n2; n3]
ind_active_views = find(~active_imgviews(1,:) & active_images);
for kk = ind_active_views,
    jj = (kk-1)*n_cam;
    active_view = find(active_imgviews(2:end, kk)')+1;
    for pp = active_view,
        kth = jj+pp;
        Rw1t = rodrigues(-Omcc(:,pp));
        Tw1 = Rw1t*(Tw_mat(:,kth)-Tcc(:,pp));
        omw1 = rodrigues(Rw1t*rodrigues(Omw_mat(:, kth)));
        if handcc(pp)~=1,
            omw1(1:2) = -omw1(1:2);
            Tw1(3) = -Tw1(3);
        end;
        Omcw_cell{kk} = [Omcw_cell{kk}, omw1];
        Tcw_cell{kk} = [Tcw_cell{kk}, Tw1];
    end;
    omw1 = Omcw_cell{kk};
    Tw1 = Tcw_cell{kk};
    n = size(omw1,2);
    switch n,
        case 0,
            fprintf(1,'All camera views of image %d are inactive!\nSet the image inactive!\n',kk);
            active_imgviews(:,kk) = 0;
            active_images(kk) = 0;
            continue;
        case 1,
            Omcw(:, kk) = omw1;
            Tcw(:, kk) = Tw1;
        otherwise,
            % transform axis angle to quaternion to compute mean rotation
            Qck = trans_quat_axis(omw1);
            Qc1 = Qck(:,1);
            for i = 2:n,
                Qcki = Qck(:,i);
                if sum(Qcki.*Qc1)<0,
                    Qck(:,i) = -Qcki;
                end;
            end;
            Q = quatmean(Qck,ones(1,n)/n);        % mean of quaternions
            if sum(Q.*Qc1)<0,
                Q = -Q;
            end;
            Tcw(:, kk) = mean(Tw1,2);        % mean of translation
            Omcw(:, kk) = trans_quat_axis(Q);
            % relative deviation of orientation and position of cameras
            deviation = max(max(sqrt(sum((Qck-Q(:,ones(1,n))).^2)),[],2), max(sqrt(sum((Tw1-Tcw(:,kk)*ones(1,n)).^2)),[],2)/norm(Tcw(:,kk)));
            if deviation>bigeps,
                fprintf(1,'The deviation value is larger than %1.0e, please stop and check variables ''Qck'' and ''Tw1''.\n', bigeps);
                flag = input('Abort calibration or not? ([]=yes, other=no) ','s');
                if isempty(flag),
                    fprintf(1,'Program terminated!\n');
                    return;
                end;
            end;
            % max deviation of camera 1's (Q,t) wrt world frame
            extrinsic_deviation(1) = max(extrinsic_deviation(1),deviation);
    end;
end;
ind_active = find(active_images);

for count=1:2,
    fprintf(1,'\nRefine extrinsic parameters: %2d\n', count);
    optim_multicams_extrinsic;
end;

flag = input('Further refine intrinsic parameters or not? ([]=no, other=yes)','s');
if isempty(flag)
    disp('Done.');
    return;
end;

for count=1:3,
    fprintf(1,'\nRefine intrinsic and extrinsic parameters: %2d\n', count);
     optim_multicams;
end;
flag = input('Change some estimation settings to further refine parameters? ([]=no, other=yes)','s');
if ~isempty(flag),
    fprintf(1,'\nCarefully choose to change estimation settings, especially distortion settings!\n');
    for pp = 1:n_cam,
        if ~isequal(est_fc_mat(:,pp),[1;1]),
            fprintf(1,'The focal length of camera %d is not fully estimized.\n',pp);
            flag = input('Do you want to change it? ([]=no, other=yes)','s');
            if ~isempty(flag),
                est_fc_mat(:,pp) = [1;1];
            end;
        end;
    end;

    if ~isequal(est_aspect_ratio_vec, ones(1,n_cam)),
        fprintf(1,'The pixel aspect ratio of all cameras are not fully estimated.\n');
        flag = input('Do you want to change it? ([]=no, other=yes)','s');
        flag = ~isempty(flag);
        while flag,
            est_aspect_ratio_vec = input(['est_aspect_ratio_vec = ([] = [' num2str(ones(1,n_cam)) '])']);
            if isempty(est_aspect_ratio_vec),
                est_aspect_ratio_vec = ones(1,n_cam);
                flag = 0;
            else
                est_aspect_ratio_vec = est_aspect_ratio_vec(:)';
                flag = length(est_aspect_ratio_vec)~=n_cam;
                if flag,
                    fprintf(1,'\nUnexpected input, please input again!\n');
                end;
            end;
        end;
    end;

    if ~isequal(center_optim_vec, ones(1,n_cam)),
        fprintf(1,'The principal point of all cameras are not fully estimated.\n');
        flag = input('Do you want to change it? ([]=no, other=yes)','s');
        flag = ~isempty(flag);
        while flag,
            center_optim_vec = input(['center_optim_vec = ([] = [' num2str(ones(1,n_cam)) '])']);
            if isempty(center_optim_vec),
                center_optim_vec = ones(1,n_cam);
                flag = 0;
            else
                center_optim_vec = center_optim_vec(:)';
                flag = length(center_optim_vec)~=n_cam;
                if flag,
                    fprintf(1,'\nUnexpected input, please input again!\n');
                end;
            end;
        end;
    end;

    if ~isequal(est_alpha_vec, ones(1,n_cam)),
        fprintf(1,'The pixel skew of all cameras are not fully estimated.\n');
        flag = input('Do you want to change it? ([]=no, other=yes)','s');
        flag = ~isempty(flag);
        while flag,
            est_alpha_vec = input(['est_alpha_vec = ([] = [' num2str(ones(1,n_cam)) '])']);
            if isempty(est_alpha_vec),
                est_alpha_vec = ones(1,n_cam);
                flag = 0;
            else
                est_alpha_vec = est_alpha_vec(:)';
                flag = length(est_alpha_vec)~=n_cam;
                if flag,
                    fprintf(1,'\nUnexpected input, please input again!\n');
                end;
            end;
        end;
    end;

    for pp = 1:n_cam,
        est_dist = est_dist_mat(:,pp);
        if ~norm(est_dist),
            fprintf(1,'\nThe distortion coefficients of camera %d are not estimized.\n',pp);
            flag = input('Do you want to estimate it? ([]=no, other=yes with cautiousness)','s');
            if ~isempty(flag),
                est_dist = input(['est_dist [k1; k2; k3; k4; k5] of camera ' num2str(pp) ': ([] = [1;1;1;1;0])])']);
                if isempty(est_dist),
                    est_dist = [1;1;1;1;0];  % by default do not estimate the 6th radial distortion
                end;
                est_dist = est_dist(:);
            end;
            n = length(est_dist);
            if n<5,
                est_dist = [est_dist;zeros(5-n,1)];
            elseif n>5,
                est_dist = est_dist(1:5);
            end;
            fprintf(1,['\nest_dist of camera ' num2str(pp) ': [%d; %d; %d; %d; %d].\n'], est_dist);
            est_dist_mat(:,pp) = est_dist;
        end;
    end;

    for count=1:2,
        fprintf(1,'\nFurther refine intrinsic and extrinsic parameters: %2d\n', count);
        optim_multicams;
    end;
    fprintf(1,'\nRefine extrinsic parameters in the end:\n');
    optim_multicams_extrinsic;
end;
disp('Done.');
