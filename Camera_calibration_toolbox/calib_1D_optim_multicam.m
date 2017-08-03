% calib_1D_optim_multicam
%
% see calib_1D_optim_mirror.
%
% Main calibration function. Computes the intrinsic and extrinsic parameters.
% Runs as a script.
%
% INPUT: x_cell: Feature locations on the images;
%        lamda: the 1D coordinates of points on the rod;
%
% OUTPUT: fc_mat: focal length of cameras;
%        cc_mat: Principal point coordinates;
%        alpha_vec: Skew coefficient;
%        kc_mat: Distortion coefficients;
%        KK: The intrinsic matrix of a camera;
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
    if exist('multicam_calib_data.mat','file')==2,
        load('multicam_calib_data.mat');
    elseif exist('multicam_simu_data.mat','file')==2,
        fprintf(1,'\nSimulation data detected! Do you want to load it?\n');
        flag = input('Load the simulation data or not? ([]=yes, other=no) ','s');
        if isempty(flag),
            load('multicam_simu_data.mat');
            clear Xrod Xori Xdir fc_mat cc_mat kc_mat alpha_vec Omcw Tcw;
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
active_images = any(active_imgviews,1);
ind_active = find(active_images);
ind_active_views = find(active_imgviews(:)');

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
        init_intrinsic_1D_multicam;
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
[~, id_mcam] = min(sum(costs, 2));
path_mcam = paths(id_mcam,:);
return;

%%% Initialization of the camera parameters for global minimization:
%%% Computes pairwise extrinsic parameters for all cameras

n_view = n_ima * n_cam;
Omw_mat = NaN(3, n_view);
Tw_mat = NaN(3, n_view);
for pp = 1:n_cam,
    fc = fc_mat(:,pp);
    cc = cc_mat(:,pp);
    kc = kc_mat(:,pp);
    alpha_c = alpha_vec(pp);
    handkk = hand_list(pp);
    active_view = find(active_imgviews(pp,:));
    for kk = active_view,
        kth = (kk-1)*n_cam+pp;
        x_kk = x_cell{kth};
        X_kk = X_cell{kth};
        if isempty(x_kk) || isnan(x_kk(1)),   % if x_kk is [] or NaN, then the image is set inactive;
            fprintf(1,'Warning: Cannot calibrate (camera %d, image %d). This image is now set inactive.\n',pp,kk);
            active_imgviews(pp,kk) = 0;
            x_cell{kth} = [];
            X_cell{kth} = [];
            dXY_mat(:, kth) = NaN(2,1);
            n_sq_mat(:, kth) = NaN(2,1);
        else
            [omwkk,Twkk] = compute_extrinsic_init2(x_kk,X_kk,fc,cc,kc,alpha_c,handkk);
            [omwkk,Twkk,JJ_kk] = compute_extrinsic_lm_refine2(x_kk,X_kk,omwkk,Twkk,handkk,fc,cc,kc,alpha_c,MaxIter2,thresh_cond);
            if check_cond && (cond(JJ_kk)> thresh_cond),
                active_imgviews(pp,kk) = 0;
                fprintf(1,'\nWarning: (camera %d, image %d) ill-conditioned. This image is now set inactive.\n',pp,kk);
            elseif any(isnan(omwkk)),
                active_imgviews(pp,kk) = 0;
                fprintf(1,'\nWarning: The extrinsic rotation is NaN! Deactivating (camera %d, image %d).\n',pp,kk);
            else
                Omw_mat(:,  kth) = omwkk;
                Tw_mat(:,  kth) = Twkk;
            end;
        end;
    end;
end;

kk = find(active_images ~= any(active_imgviews,1));
if ~isempty(kk),
    desactivated_images = [desactivated_images kk];
    fprintf(1,'WARNING: Cannot calibrate all cameras of image %d.\n',kk);
    fprintf(1,'Set active_images(%d)=0;\n',kk);
end;
active_images = any(active_imgviews,1);

% keyboard;

%% ------------------------------------------ Main Optimization:

fprintf(1,'\nMain calibration optimization procedure - Bundle adjustment with %d cameras\n', n_cam);
fprintf(1,'Sparse Levenberg-Marquardt iterations: ');

%%% Initialization of the global parameter vector:
ncam10 = 10*n_cam;
nview6 = 6*n_view;
intrinsic_param = reshape([fc_mat; cc_mat; alpha_vec; kc_mat],ncam10,1);
extrinsic_param = reshape([Omw_mat; Tw_mat],nview6,1);

intr_update = intrinsic_param;
extr_update = extrinsic_param;
ex = []; % Global error vector
for pp=1:n_cam,
    fc = fc_mat(:,pp);
    cc = cc_mat(:,pp);
    kc = kc_mat(:,pp);
    alpha_c = alpha_vec(pp);
    handkk = hand_list(pp);
    active_view = find(active_imgviews(pp,:));
    for kk = active_view,
        % Reproject the patterns on the images, and compute the pixel errors:
        kth = (kk-1)*n_cam+pp;
        omwkk = Omw_mat(:, kth);
        Twkk = Tw_mat(:, kth);
        X_kk = X_cell{kth};
        x_kk = x_cell{kth};
        x = project_points_mirror2(X_kk,omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
        ex_kk = x_kk - x;
        ex = [ex, ex_kk];
    end;
end;
err_std0 = std(ex,0,2);

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
        U = sparse([],[],[],ncam10,ncam10,10*ncam10);
        V = sparse([],[],[],nview6,nview6,6*nview6);
        W = sparse([],[],[],ncam10,nview6,ncam10*nview6);

        ea = zeros(ncam10,1);        % A'*ex
        eb = zeros(nview6,1);         % B'*ex
        for pp=1:n_cam,
            ii = (pp-1)*10;
            handkk = hand_list(pp);
            fc = intrinsic_param(ii+1 : ii+2);
            cc = intrinsic_param(ii+3 : ii+4);
            alpha_c = intrinsic_param(ii+5);
            kc = intrinsic_param(ii+6 : ii+10);
            est_aspect_ratio = est_aspect_ratio_vec(pp);
            ind_active_views = find(active_imgviews(pp,:));
            for kk = ind_active_views,
                kth = (kk-1)*n_cam+pp;
                jj = (kth-1)*6;
                omwkk = extrinsic_param(jj+1 : jj+3);
                Twkk = extrinsic_param(jj+4 : jj+6);
                if any(isnan(omwkk)),
                    fprintf(1,'Extrinsic parameters of (camera %d, image %d) do not exist!\n',pp,kk);
                    return;
                end;
                X_kk = X_cell{kth};
                x_kk = x_cell{kth};
                if est_aspect_ratio,
                    [x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points_mirror2(X_kk,omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
                else
                    [x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points_mirror2(X_kk,omwkk,Twkk,handkk,fc(1),cc,kc,alpha_c);
                    dxdf = repmat(dxdf,[1 2]);
                end;
                ex_kk = x_kk - x;
                Akk = [dxdf dxdc dxdalpha dxdk];
                Bkk = [dxdom dxdT];

                U(ii+1 : ii+10, ii+1 : ii+10) = U(ii+1 : ii+10, ii+1 : ii+10) + sparse(Akk'*Akk);
                ea(ii+1 : ii+10) = ea(ii+1 : ii+10) + Akk'*ex_kk(:);
                % Check if this view is ill-conditioned:
                if check_cond && (cond(Bkk)>thresh_cond),
                    fprintf(1,['\nWarning: (camera %d, image %d) ill-conditioned. This view is now set inactive. \n' ...
                                   '(note: to disactivate this option, set check_cond=0)\n'],pp,kk);
                    active_imgviews(pp,kk) = 0;
                    extrinsic_param(jj+1 : jj+6) = NaN(6,1);
                else
                    V(jj+1 : jj+6, jj+1 : jj+6) = sparse(Bkk'*Bkk);
                    W(ii+1 : ii+10, jj+1 : jj+6) = sparse(Akk'*Bkk);
                    eb(jj+1 : jj+6) = Bkk'*ex_kk(:);
                end;
            end;
        end;

        active_view = active_imgviews(:)';
        % The following vector helps to select the variables to update (for only active images):
        selected_invar = [est_fc_mat; ones(2,1)*center_optim_vec; est_alpha_vec; est_dist_mat];
        selected_invar(2,:) = selected_invar(2,:).*(est_aspect_ratio_vec | ~est_fc_mat(1,:));
        selected_invar = selected_invar(:)';
        selected_exvar = reshape(ones(6,1)*active_view,1,nview6);

        ind_va = find(selected_invar);
        ind_vb = find(selected_exvar);
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
    extr_innov = V_lm\(eb-W'*intr_innov);             % db

    intr_update(ind_va) = intrinsic_param(ind_va) + intr_innov;     % updated parameters
    extr_update(ind_vb) = extrinsic_param(ind_vb) + extr_innov;

    % compute reprojection error vector
    ex = [];
    for pp=1:n_cam,
        ii = (pp-1)*10;
        handkk = hand_list(pp);
        % New intrinsic parameters:
        if ~est_aspect_ratio_vec(pp) && all(est_fc_mat(:,pp)),
            intr_update(ii+2) = intr_update(ii+1);
        end;
        fc = intr_update(ii+1 : ii+2);
        cc = intr_update(ii+3 : ii+4);
        alpha_c = intr_update(ii+5);
        kc = intr_update(ii+6 : ii+10);
        if center_optim_vec(pp) && (cc(1)<-.5 || cc(1)>imsize(1,pp)-.5 || cc(2)<-.5 || cc(2)>imsize(2,pp)-.5),
            fprintf(1,'Warning: it appears that the principal point of camera %d cannot be estimated.\n', pp);
            center_optim_vec(pp) = 0;
            cc = cc_mat(:,pp);
            intr_update(ii+3 : ii+4) = cc;
        end;
        ind_active_views = find(active_imgviews(pp,:));
        for kk = ind_active_views,
            kth = (kk-1)*n_cam+pp;
            jj = (kth-1)*6;
            omwkk = extr_update(jj+1 : jj+3);
            Twkk = extr_update(jj+4 : jj+6);
            X_kk = X_cell{kth};
            x_kk = x_cell{kth};
            x = project_points_mirror2(X_kk,omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
            ex_kk = x_kk - x;
            ex = [ex, ex_kk];
        end;
    end;
    ex = ex(:);
    ex2_update = dot(ex,ex);
    if ex2_update < ex2,
        lamda = lamda/10;
        intrinsic_param = intr_update;
        extrinsic_param = extr_update;
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

kk = find(active_images ~= any(active_imgviews,1));
if ~isempty(kk),
    desactivated_images = [desactivated_images kk];
    fprintf(1,'WARNING: Cannot calibrate all camera views of image %d.\n',kk);
    fprintf(1,'Set active_images(%d)=0;\n',kk);
end;
if ~isempty(desactivated_images),
    fprintf(1,['List of images left desactivated: ' num2str(desactivated_images) '.\n' ] );
end;
active_images = any(active_imgviews,1);
ind_active = find(active_images);

%%%--------------------------- Computation of the error of estimation:

fprintf(1,'\nEstimation of uncertainties...');

% Extraction of the paramters for computing the reprojection error:
intrinsic_param = reshape(intrinsic_param, 10, n_cam);
fc_mat = intrinsic_param(1:2,:);
cc_mat = intrinsic_param(3:4,:);
alpha_vec = intrinsic_param(5,:);
kc_mat = intrinsic_param(6:10,:);

extrinsic_param = reshape(extrinsic_param, 6, n_view);
Omw_mat = extrinsic_param(1:3,:);
Tw_mat = extrinsic_param(4:6,:);

% Reproject the patterns on the images, and compute the pixel errors:
y_cell = cell(1, n_view);  % Reprojected points
ex_cell = y_cell;             % Reprojected error
H_cell = y_cell;
err_cam = zeros(2,n_cam);
ex = [];
for pp = 1:n_cam,
    fc = fc_mat(:,pp);
    cc = cc_mat(:,pp);
    alpha_c = alpha_vec(pp);
    kc = kc_mat(:,pp);
    handkk = hand_list(pp);
    % Calibration matrix:
    KK = [fc(1) fc(1)*alpha_c cc(1);0 fc(2) cc(2); 0 0 1];
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
        Rwkk = rodrigues(omwkk);
        Hkk = KK * [Rwkk(:,1) Rwkk(:,2) Twkk];
        H_cell{kth} = Hkk / Hkk(3,3);
    end;
    err_cam(:,pp) = std(cell2mat(ex_cell(pp:n_cam:end)),0,2);
end;

err_std = std(ex,0,2);
sigma_x = std(ex(:));

% Compute the the standard deviation of parameters from std(ex):
% ex = X-f(P),  cov(param,param) = inv((JJ'* inv(cov(ex,ex))* JJ))
JJ2 = [U, W; W', V];    % JJ2 = JJ'*JJ
JJ2_inv = speye(size(JJ2))/JJ2;

% Extraction of the final intrinsic paramaters:
param_error = zeros(ncam10,1);           % take value as perfect if not optimized (no error).
nva = length(ind_va);
param_error(ind_va) = 3*sqrt(abs(full(diag(JJ2_inv(1:nva,1:nva)))))*sigma_x;     % 3 sigma principle
param_error = reshape(param_error,10,n_cam);
ind_valid = find(~est_aspect_ratio_vec & all(est_fc_mat,1));
param_error(2,ind_valid) = param_error(1,ind_valid);

fc_mat_error = param_error(1:2,:);
cc_mat_error = param_error(3:4,:);
alpha_vec_error = param_error(5,:);
kc_mat_error = param_error(6:10,:);

% Extraction of the final extrinsic paramaters:
param_error = zeros(nview6,1);           % take value as perfect if not optimized (no error).
nvb = length(ind_vb);
param_error(ind_vb) = 3*sqrt(abs(full(diag(JJ2_inv(nva+1:nva+nvb, nva+1:nva+nvb)))))*sigma_x;     % 3 sigma principle
param_error = reshape(param_error,6,n_view);
Omw_mat_error = param_error(1:3,:);
Tw_mat_error = param_error(4:6,:);

fprintf(1,'Done with the uncertainties!\n');

for pp = 1:n_cam,
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
    Omcw(:, kk) = omw1;
    Tcw(:, kk) = Tw1;
    active_view = find(active_imgviews(2:end, kk)')+1;
    for pp = active_view,
        kth = jj+pp-1;
        omwkk = Omw_mat(:, kth);
        Twkk = Tw_mat(:, kth);
        handkk = handcc(pp);
        Rw1t = rodrigues(-omw1);
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
