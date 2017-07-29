%go_calib_optim_mirror
%
% version with pure LM algorithm
% see go_calib_optim_mirror2 using axis angle rotation
%
%Main calibration function. Computes the intrinsic and extrinsic parameters.
%Runs as a script.
%
%INPUT:x_cell: Feature locations on the images
%       X_cell: Corresponding grid coordinates
%
%OUTPUT: fc: Camera focal length
%        cc: Principal point coordinates
%        alpha_c: Skew coefficient
%        kc: Distortion coefficients
%        KK: The camera matrix (containing fc and cc)
%        Qw_mat: quaternion rotation vectors attached to the grid positions in space
%        Tw_mat: 3D translation vectors attached to the grid positions in space
%
%Method: Minimizes the pixel reprojection error in the least squares sense over the intrinsic
%        camera parameters, and the extrinsic parameters (3D locations of the grids in space)
%
%Note: If the intrinsic camera parameters (fc, cc, kc) do not exist before, they are initialized through
%      the function init_intrinsic_mirror.m. Otherwise, the variables in memory are used as initial guesses.
%
%Note: The row vector active_images consists of zeros and ones. To deactivate an image, set the
%      corresponding entry in the active_images vector to zero.
%
%VERY IMPORTANT: This function works for 2D and 3D calibration rigs, except for init_intrinsic_mirror.m
%that is so far implemented to work only with 2D rigs.
%In the future, a more general function will be there.
%For now, if using a 3D calibration rig, focal length will be initialized with 35 degree of FOV angle.

if ~exist('x_cell','var'),
    if exist('multimirror_calib_data.mat','file')==2,
        load('multimirror_calib_data.mat');
    else
        fprintf(1,['\nThere is no data required for calibration! \nPlease run "multimirrors_gui" '...
            'and click the first two buttons!\n"Read images" , then "Extract grid corners"...\n']);
        return;
    end;
end;

check_extracted_mirror;
ind_active_views = find(active_imgviews(:)');

if ~exist('desactivated_images','var'),
    desactivated_images = [];
end;

if ~exist('est_aspect_ratio','var'),
    fprintf(1,'Do you want to estimate the aspect ratio of pixel?\n Set to zero if you don''t!\n');
    est_aspect_ratio = input('est_aspect_ratio = ([] = 1)');
    if isempty(est_aspect_ratio),
        est_aspect_ratio = 1;
    end;
end;

if ~exist('est_fc','var'),
    fprintf(1,'Do you want to estimate the focal length?\n Set to zero if you don''t!\n');
    est_fc = input('est_fc [fc1; fc2] = ([] = [1;1])');
    if isempty(est_fc),
        est_fc = [1;1];
    end;
end;

if ~exist('MaxIter','var'),
    MaxIter = 30; % Maximum number of iterations in the main LM algorithm
end;

if ~exist('MaxIter2','var'),
    MaxIter2 = 20; % Maximum number of iterations in the function compute_extrinsic_lm_refine
end;

if ~exist('check_cond','var'),
    check_cond = 0; % Set this variable to 0 in case you don't want to extract view dynamically
end;

if ~exist('center_optim','var'),
    fprintf(1,'Do you want to estimate the principal point?\n Set to zero if you don''t!\n');
    center_optim = input('center_optim = ([] = 1)');
    if isempty(center_optim),
        center_optim = 1;
    end;
end;

if ~exist('est_dist','var'),
    fprintf(1,'Do you want to estimate the distortion coefficients?\n Set to zero if you don''t!\n');
    est_dist = input('est_dist [k1; k2; k3; k4; k5] = ([] = [1;1;1;1;0])');
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


if ~exist('est_alpha','var'),
    fprintf(1,'Do you want to estimate the skew of pixels?\n Set to zero if you don''t!\n');
    est_alpha = input('est_alpha = ([] = 0)');
    if isempty(est_alpha),
        est_alpha = 0; % by default do not estimate skew
    end;
end;

% Little fix in case of stupid values in the binary variables:
center_optim = double(~~center_optim);
est_alpha = double(~~est_alpha);
est_dist = double(~~est_dist);
est_fc = double(~~est_fc);
est_aspect_ratio = double(~~est_aspect_ratio);


fprintf(1,'\n');

if ~exist('nx','var') || ~exist('ny','var'),
    fprintf(1,'WARNING: No image size (nx,ny) available. Please set the image size manually!\n');
    nx = input('nx = ([] = 640)');
    if isempty(nx),
        nx = 640;
    end;
    ny = input('ny = ([] = 480)');
    if isempty(ny),
        ny = 480;
    end;
end;

% Check 3D-ness of the calibration rig:
rig3D = 0;
for kk = ind_active_views,
    X_kk = X_cell{kk};
    if is3D(X_kk),
        rig3D = 1;
    end;
end;

if center_optim && (length(ind_active_views) < 2) && ~rig3D,
    fprintf(1,'WARNING: Principal point rejected from the optimization when using one image and planar rig (center_optim = 1).\n');
    center_optim = 0; %%% when using a single image, please, no principal point estimation!!!
    est_alpha = 0;
end;

if center_optim && (length(ind_active_views) < 5) && ~rig3D,
    fprintf(1,'WARNING: The principal point estimation may be unreliable (using less than 5 images for calibration).\n\n');
end;

if ~est_aspect_ratio,
    fprintf(1,'Aspect ratio not optimized (est_aspect_ratio = 0) -> fc(1)=fc(2). Set est_aspect_ratio to 1 for estimating aspect ratio.\n\n');
else
    if isequal(est_fc,[1;1]),
        fprintf(1,'Aspect ratio optimized (est_aspect_ratio = 1) -> both components of fc are estimated (DEFAULT).\n\n');
    end;
end;

if ~isequal(est_fc,[1;1]),
    if isequal(est_fc,[1;0]),
        fprintf(1,'The first component of focal (fc(1)) is estimated, but not the second one (est_fc=[1;0])\n\n');
    else
        if isequal(est_fc,[0;1]),
            fprintf(1,'The second component of focal (fc(1)) is estimated, but not the first one (est_fc=[0;1])\n\n');
        else
            fprintf(1,'The focal vector fc is not optimized (est_fc=[0;0])\n\n');
        end;
    end;
end;

if ~prod(double(est_dist)),
    fprintf(1,'Distortion not fully estimated (defined by the variable est_dist):\n');
    if ~est_dist(1),
        fprintf(1,'     Second order distortion not estimated (est_dist(1)=0).\n');
    end;
    if ~est_dist(2),
        fprintf(1,'     Fourth order distortion not estimated (est_dist(2)=0).\n');
    end;
    if ~est_dist(5),
        fprintf(1,'     Sixth order distortion not estimated (est_dist(5)=0) - (DEFAULT) .\n');
    end;
    if ~prod(double(est_dist(3:4))),
        fprintf(1,'     Tangential distortion not estimated (est_dist(3:4)~=[1;1]).\n');
    end;
end;

if ~exist('fc','var') && rig3D,
    FOV_angle = 35; %field of view in degrees: for135 camera, 35 degree of FOV is about 57 mm focal length。
    fprintf(1,['Initialization of the focal length to a FOV of ' num2str(FOV_angle) ' degrees.\n\n']);
    fc = (nx/2)/tan(pi*FOV_angle/360) * ones(2,1);    % FOV_angle=2*atan(nx/(2*fc))
end;

if ~exist('fc','var'),
    % Initialization of the intrinsic parameters with no distortion:
    fprintf(1,'Initialization of the intrinsic parameters using the constraints of homography.\n\n');
    init_intrinsic_mirror;    % see also init_intrinsic_mirror2;
end;

% Initialization of the intrinsic parameters (if necessary: rig3D==1)
if ~exist('cc','var'),
    fprintf(1,'\nInitialization of the principal point at the center of the image.\n');
    cc = [(nx-1)/2;(ny-1)/2];
end;

if ~exist('kc','var'),
    fprintf(1,'Initialization of the image distortion to zero.\n');
    kc = zeros(5,1);
else
    n = length(kc);
    if n<5,
        fprintf(1,'\nRadial distortion model up to the 6th degree.\n\n');
        kc = [kc;zeros(5-n,1)];
    elseif n>5,
        kc = kc(1:5);
    end;
end;

if ~center_optim, % In the case where the principal point is not estimated, keep it at the center of the image
    fprintf(1,'Principal point not optimized (center_optim=0). ');
    if exist('cc','var'),
        fprintf(1,'Note: to set it in the middle of the image, clear variable cc, and run calibration again.\n\n');
    end;
else
    fprintf(1,'Principal point optimized (center_optim=1) - (DEFAULT). To reject principal point, set center_optim=0\n\n');
end;

if ~center_optim && est_alpha,
    fprintf(1,'WARNING: Since there is no principal point estimation (center_optim=0), no skew estimation (est_alpha = 0)\n\n');
    est_alpha = 0;
end;

if ~est_alpha,
    fprintf(1,'Skew not optimized (est_alpha=0) - (DEFAULT)\n\n');
    alpha_c = 0;
else
    fprintf(1,'Skew optimized (est_alpha=1). To disable skew estimation, set est_alpha=0.\n\n');
end;

if ~exist('alpha_c','var'),
    fprintf(1,'Initialization of the image skew to zero.\n');
    alpha_c = 0;
end;

if ~prod(double(est_dist)),
    % If no distortion estimated, set to zero the variables that are not estimated
    kc = kc .* est_dist;
end;

if ~prod(double(est_fc)),
    fprintf(1,'Warning: The focal length is not fully estimated (est_fc ~= [1;1])\n');
end;

% Conditioning threshold for view rejection
thresh_cond = 1e8;
% threshold to terminate the main LM iteration
gradeps = eps*1e8;

%%% Initialization of the extrinsic parameters for global minimization:
%%% Computes the extrinsic parameters for all the active calibration images

n_view = n_ima * n_cam;
Qw_mat = NaN(4, n_view);
Tw_mat = NaN(3, n_view);
for pp = 1:n_cam,
    handkk = hand_list(pp);
    active_view = find(active_imgviews(pp,:));
    for kk = active_view,
        kth = (kk-1)*n_cam+pp;
        x_kk = x_cell{kth};
        X_kk = X_cell{kth};
        if isempty(x_kk) || isnan(x_kk(1)),   % 若x_kk为[]或NaN，则此视角无效;
            fprintf(1,'Warning: Cannot calibrate with view %d of image %d. This view is now set inactive.\n',pp,kk);
            active_imgviews(pp,kk) = 0;
            x_cell{kth} = [];
            X_cell{kth} = [];
            dXY_mat(:, kth) = NaN(2,1);
            n_sq_mat(:, kth) = NaN(2,1);
        else
            [Qwkk,Twkk] = compute_extrinsic_init_mirror(x_kk,X_kk,fc,cc,kc,alpha_c,handkk);
            [Qwkk,Twkk,JJ_kk] = compute_extrinsic_lm_refine(x_kk,X_kk,Qwkk,Twkk,handkk,fc,cc,kc,alpha_c,MaxIter2,thresh_cond);
            if check_cond && (cond(JJ_kk)> thresh_cond),
                active_imgviews(pp,kk) = 0;
                fprintf(1,'\nWarning: (view %d, image %d) ill-conditioned. The view is now set inactive.\n',pp,kk);
            elseif any(isnan(Qwkk)),
                active_imgviews(pp,kk) = 0;
                fprintf(1,'\nWarning: The extrinsic rotation is NaN! Deactivating (view %d, image %d).\n',pp,kk);
            else
                Qw_mat(:,  kth) = Qwkk;
                Tw_mat(:,  kth) = Twkk;
            end;
        end;
    end;
end;
kk = find(active_images ~= any(active_imgviews,1));
if ~isempty(kk),
    desactivated_images = [desactivated_images kk];
    fprintf(1,'WARNING: Cannot calibrate all views of image %d.\n',kk);
    fprintf(1,'Set active_images(%d)=0;\n',kk);
end;
active_images = any(active_imgviews,1);
ind_active_views = find(active_imgviews(:)');

% keyboard;

%------------------------------------------ Main Optimization:

fprintf(1,'\nMain calibration optimization procedure - Number of views: %d\n',length(ind_active_views));
fprintf(1,'Sparse Levenberg-Marquardt iterations: ');

%%% Initialization of the global parameter vector:
nview7 = 7*n_view;
intrinsic_param = [fc;cc;alpha_c;kc];
extrinsic_param = reshape([Qw_mat; Tw_mat],nview7,1);

intr_update = intrinsic_param;
extr_update = extrinsic_param;
ex = []; % Global error vector
for kk = ind_active_views,
    omwkk = Omw_mat(:, kk);
    Twkk = Tw_mat(:, kk);
    % Reproject the patterns on the images, and compute the pixel errors:
    X_kk = X_cell{kk};
    x_kk = x_cell{kk};
    handkk = hand_list(mod(kk-1,n_cam)+1);
    x = project_points_mirror2(X_kk,omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
    ex_kk = x_kk - x;
    ex = [ex, ex_kk];
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
        U = sparse([],[],[],10,10,100);
        V = sparse([],[],[],nview7,nview7,7*nview7);
        W = sparse([],[],[],10,nview7,10*nview7);
        ea = zeros(10,1);            % A'*ex
        eb = zeros(nview7,1);      % B'*ex

        f = intrinsic_param(1:2);
        c = intrinsic_param(3:4);
        alpha = intrinsic_param(5);
        k = intrinsic_param(6:10);

        for kk = ind_active_views,
            ii = 7*(kk-1);
            Qwkk = extrinsic_param(ii+1 : ii+4);
            Twkk = extrinsic_param(ii+5 : ii+7);
            if any(isnan(Qwkk)),
                fprintf(1,'Extrinsic parameters at view %d of frame %d do not exist!\n',mod(kk-1,n_cam)+1,ceil(kk/n_cam));
                return;
            end;
            X_kk = X_cell{kk};
            x_kk = x_cell{kk};
            handkk = hand_list(mod(kk-1,n_cam)+1);
            if est_aspect_ratio,
                [x,dxdQ,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points_mirror(X_kk,Qwkk,Twkk,handkk,f,c,k,alpha);
            else
                [x,dxdQ,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points_mirror(X_kk,Qwkk,Twkk,handkk,f(1),c,k,alpha);
                dxdf = repmat(dxdf,[1 2]);
            end;
            ex_kk = x_kk - x;
            Akk = [dxdf dxdc dxdalpha dxdk];
            Bkk = [dxdQ dxdT];

            U = U + sparse(Akk'*Akk);
            ea = ea + Akk'*ex_kk(:);
            % Check if this view is ill-conditioned:
            if check_cond && (cond(Bkk)> thresh_cond),
                fprintf(1,['\nWarning: View %d of frame %d ill-conditioned. This view is now set inactive. \n' ...
                    '(note: to disactivate this option, set check_cond=0)\n'],mod(kk-1,n_cam)+1,ceil(kk/n_cam));
                active_imgviews(kk) = 0;
                extrinsic_param(ii+1 : ii+7) = NaN(7,1);
            else
                V(ii+1 : ii+7, ii+1 : ii+7) = sparse(Bkk'*Bkk);
                W(1:10, ii+1 : ii+7) = sparse(Akk'*Bkk);
                eb(ii+1 : ii+7) = Bkk'*ex_kk(:);
            end;
        end;
        active_view = active_imgviews(:)';
        ind_active_views = find(active_view);

        % The following vector helps to select the variables to update (for only active images):
        selected_invar = [est_fc; center_optim*ones(2,1); est_alpha; est_dist]';
        selected_exvar = reshape(ones(7,1)*active_view,1,nview7);
        if ~est_aspect_ratio && est_fc(1),
            selected_invar(2) = 0;
        end;

        ind_va = find(selected_invar);
        ind_vb = find(selected_exvar);

        U = U(ind_va,ind_va);
        V = V(ind_vb,ind_vb);
        W = W(ind_va,ind_vb);

        ea = ea(ind_va);
        eb = eb(ind_vb);
    end;

    U_lm = U + diag(lamda*diag(U));     % U + lamda*speye(size(U));
    V_lm =V + diag(lamda*diag(V));      % V + lamda*speye(size(V));
    Y = W/V_lm;

    intr_innov = (U_lm-Y*W')\(ea-Y*eb);              % da
    extr_innov = V_lm\(eb-W'*intr_innov);             % db
    intr_update(ind_va) = intrinsic_param(ind_va) + intr_innov;     % updated parameters
    extr_update(ind_vb) = extrinsic_param(ind_vb) + extr_innov;

    % New intrinsic parameters:
    if ~est_aspect_ratio && all(est_fc),
        intr_update(2) = intr_update(1);
    end;
    f = intr_update(1:2);
    c = intr_update(3:4);
    alpha = intr_update(5);
    k = intr_update(6:10);

    if center_optim && (c(1)<0 || c(1)>nx || c(2)<0 || c(2)>ny),
        fprintf(1,'Warning: it appears that the principal point cannot be estimated. Setting center_optim = 0\n');
        center_optim = 0;
        c = cc;
        intr_update(3:4) = cc;
    end;

    % compute reprojection error vector
    ex = [];
    for kk = ind_active_views,
        ii = 7*(kk-1);
        Qwkk = extr_update(ii+1 : ii+4);
        Qwkk = Qwkk/norm(Qwkk);
        extr_update(ii+1 : ii+4) = Qwkk;
        Twkk = extr_update(ii+5 : ii+7);
        X_kk = X_cell{kk};
        x_kk = x_cell{kk};
        handkk = hand_list(mod(kk-1,n_cam)+1);
        x = project_points_mirror(X_kk,Qwkk,Twkk,handkk,f,c,k,alpha);
        ex_kk = x_kk - x;
        ex = [ex, ex_kk];
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
    fprintf(1,'WARNING: Cannot calibrate all views of image %d.\n',kk);
    fprintf(1,'Set active_images(%d)=0;\n',kk);
end;
if ~isempty(desactivated_images),
    fprintf(1,['List of images left desactivated: ' num2str(desactivated_images) '.\n']);
end;
active_images = any(active_imgviews,1);
ind_active = find(active_images);

%%%--------------------------- Computation of the error of estimation:

fprintf(1,'\nEstimation of uncertainties...');

% Extraction of the paramters for computing the reprojection error:
fc = intrinsic_param(1:2);
cc = intrinsic_param(3:4);
alpha_c = intrinsic_param(5);
kc = intrinsic_param(6:10);
% Calibration matrix:
KK = [fc(1) fc(1)*alpha_c cc(1);0 fc(2) cc(2); 0 0 1];

% Reproject the patterns on the images, and compute the pixel errors:
extrinsic_param = reshape(extrinsic_param, 7, n_view);
Qw_mat = extrinsic_param(1:4,:);
Tw_mat = extrinsic_param(5:7,:);

y_cell = cell(1, n_view);     % Reprojected points
ex_cell = y_cell;             % Reprojected error
H_cell = y_cell;              % recomputer the collineations
ex = [];
ind_active_views = find(active_imgviews(:)');
for kk = ind_active_views,
    Qwkk = Qw_mat(:, kk);
    Twkk = Tw_mat(:, kk);
    X_kk = X_cell{kk};
    x_kk = x_cell{kk};
    handkk = hand_list(mod(kk-1,n_cam)+1);
    y_kk = project_points_mirror(X_kk,Qwkk,Twkk,handkk,fc,cc,kc,alpha_c);
    ex_kk = x_kk- y_kk;
    ex = [ex, ex_kk];
    y_cell{kk} = y_kk;
    ex_cell{kk} = ex_kk;
    Rwkk = trans_quat_mat(Qwkk);
    Hkk = KK * [Rwkk(:,1) Rwkk(:,2) Twkk];
    H_cell{kk} = Hkk / Hkk(3,3);
end;

err_std = std(ex,0,2);
sigma_x = std(ex(:));

% Compute the the standard deviation of parameters from std(ex):
% ex = X-f(P),  cov(param,param) = inv((JJ'* inv(cov(ex,ex))* JJ))
JJ2 = [U, W; W', V]; % not bad for sparse matrices!!
JJ2_inv = speye(size(JJ2))/JJ2;
param_error = zeros(nview7+10,1);        % take value as perfect if not optimized (no error).
param_error([ind_va,ind_vb+10]) = 3*sqrt(abs(full(diag(JJ2_inv))))*sigma_x;     % 3 sigma principle
if ~est_aspect_ratio && all(est_fc),
    param_error(2) = param_error(1);
end;

%%% Extraction of the final intrinsic and extrinsic paramaters:
fc_error = param_error(1:2);
cc_error = param_error(3:4);
alpha_c_error = param_error(5);
kc_error = param_error(6:10);

% Extract the extrinsic paramters
param_error = reshape(param_error(11 : nview7+10),7,n_view);
Qw_mat_error = param_error(1:4,:);
Tw_mat_error = param_error(5:7,:);

fprintf(1,'Done with the uncertainties!\n');

err_cam = zeros(2,n_cam);
fprintf(1,'\n\nCalibration results after optimization (with uncertainties):\n\n');
fprintf(1,'Focal Length:        fc = [ %3.5f   %3.5f ] ? [ %3.5f   %3.5f ]\n',[fc;fc_error]);
fprintf(1,'Principal point:     cc = [ %3.5f   %3.5f ] ? [ %3.5f   %3.5f ]\n',[cc;cc_error]);
fprintf(1,'Skew:           alpha_c = %3.5f  ?  %3.5f  =>  angle of pixel axes = %3.5f  ?  %3.5f degrees\n', [alpha_c; alpha_c_error; 90-atan(alpha_c)*180/pi; atan(alpha_c_error)*180/pi]);
fprintf(1,'Distortion:          kc = [ %3.5f   %3.5f   %3.5f   %3.5f   %3.5f ] ? [ %3.5f   %3.5f   %3.5f   %3.5f   %3.5f ]\n',[kc;kc_error]);
fprintf(1,'Pixel error of every camera view:\n');
for pp = 1:n_cam,
    err_cam(:,pp) = std(cell2mat(ex_cell(pp:n_cam:end)),0,2);
    fprintf(1,'   view %d:         err = [ %3.5f   %3.5f ]\n',pp,err_cam(:,pp));
end;
fprintf(1,'\nNote: The pixel error are approximately three times the standard deviations (for reference).\n\n\n');

%%% Some recommendations to the user to reject some of the difficult unkowns... Still in debug mode.
alpha_c_min = alpha_c - alpha_c_error/2;
alpha_c_max = alpha_c + alpha_c_error/2;
if (alpha_c_min < 0) && (alpha_c_max > 0),
    fprintf(1,'Note: the skew coefficient alpha_c is found to be equal to zero (within its uncertainty).\n');
    fprintf(1,'Setting est_alpha=0;\n\n');
    est_alpha=0;
end;

kc_min = kc - kc_error/2;
kc_max = kc + kc_error/2;
prob_kc = (kc_min < 0) & (kc_max > 0);
if ~(prob_kc(3) && prob_kc(4))
    prob_kc(3:4) = [0;0];
end;
if any(prob_kc),
    est_dist = est_dist & ~prob_kc;
    fprintf(1,'Note: Some distortion coefficients are found equal to zero (within their uncertainties).\n');
    fprintf(1,'Setting est_dist=[%d;%d;%d;%d;%d];\n\n',est_dist);
end;
fprintf(1,['For accurate and stable error estimates, it is recommended to run Calibration twice.\n' ...
    'You may using projection to recompute the conners to improve the quality of your data!\n\n']);


%% Refine the camera parameters if more than one view is available in one image.
% First,determine the relative orientation and position of each view to the first one.
% Then use intrinsic parameters, relations between cameras, and the extrinsic parameters
% of the first view as the initial value for the refinement.

refine_multicam = 0;
if n_cam>1,
    fprintf(1,['More than one view detected, do you want to trigger the multi-view refinement now?\n'...
                   'If you have already run the Calibration twice, it is recommended to trigger refinement!\n']);
    refine_multicam = input('Refine the multi-view calibration or not? ([]=yes, other=no) ','s');
    refine_multicam = isempty(refine_multicam);
end;
if ~refine_multicam,
    return;
end;

% Initialization of the extrinsic parameters for multi-view refinement:
hand1 = hand_list(1);
handcc = hand1*hand_list;        % handcc表示视角之间的手性关系
Qcw = NaN(4, n_ima);             % Qcw和Tcw表示第一视角与场景参照物之间的方位关系
Tcw = NaN(3, n_ima);

Qcc = NaN(4, n_cam);             % Qcc和Tcc表示视角之间的方位关系(view pp wrt view 1)
Qcc(:,1) = [zeros(3,1);1];
Tcc = NaN(3, n_cam);
Tcc(:,1) = zeros(3,1);

Qcc_cell = cell(1, n_cam);       % Qcc_cell和Tcc_cell表示由不同参照物确定的视角间的方位关系（待处理）
Tcc_cell = Qcc_cell;
ind_active_views = find(active_imgviews(1,:));

% % For quaternion Q=[q1;q2;q3;q4] and Hzkk = diag([1,1, -1]),
% % if R=Hzkk*trans_quat_mat(Q)*Hzkk, let Q1= trans_quat_mat(R).
% % Then we have  Q1=[-q1;-q2;q3;q4];
for kk = ind_active_views,
    jj = (kk-1)*n_cam+1;
    Qw1 = Qw_mat(:, jj);              %   Rc1= trans_quat_mat(Qw1);
    Tw1= Tw_mat(:, jj);
    Qcw(:, kk) = Qw1;
    Tcw(:, kk) = Tw1;
    active_view = find(active_imgviews(2:end, kk)')+1;
    for pp = active_view,
        kth = jj+pp-1;
        Qwkk = Qw_mat(:, kth);               % Rwkk = trans_quat_mat(Qwkk);
        Twkk = Tw_mat(:, kth);
        handkk = handcc(pp);
        Qw1t = Qw1;
        if handkk==1,                        % Rwkk = Rwkk*Rc1';  Tp = Twkk - Rwkk*Tw1;
            Qw1t(1:3) = -Qw1t(1:3);          % invQw1=[-Qw1(1:3); Qw1(4)]
        else                                 % Rwkk = Rwkk*Hzkk*Rc1'*Hzkk;  Tp = Twkk - Rwkk*Hzkk*Tw1;
            Qw1t(3) = -Qw1t(3);              % trans_quat_mat(Hzkk*Rc1'*Hzkk)
        end;
        Qp = quatmul(Qwkk, Qw1t);
        Tp = rigid_trans(-Tw1,Qp,Twkk,handkk);
        Qcc_cell{pp} = [Qcc_cell{pp}, Qp];
        Tcc_cell{pp} = [Tcc_cell{pp}, Tp];
    end;
end;

% relative deviation checking of close rotation and translation of each view wrt the 1st view.
bigeps = 1e-2;
extrinsic_deviation = zeros(1,n_cam);
for pp=2:n_cam,
    Qp = Qcc_cell{pp};
    Tp = Tcc_cell{pp};
    n = size(Qp,2);
    switch n,
        case 0,
            fprintf(1,'No transformation of view %d wrt view 1 found!\n',pp);
            active_imgviews(pp,:) = 0;
            continue;
        case 1,
            Qcc(:, pp) = Qp;
            Tcc(:, pp) = Tp;
        otherwise,
            Qw1 = Qp(:,1);
            for i = 2:n,
                Qwkk = Qp(:,i);
                if sum(Qwkk.*Qw1)<0,
                    Qp(:,i) = -Qwkk;
                end;
            end;
            Q = quatmean(Qp,ones(1,n)/n);        % mean of quaternions
            if sum(Q.*Qw1)<0,
                Q = -Q;
            end;
            Qcc(:, pp) = Q;
            Tcc(:, pp) = mean(Tp,2);        % mean of translation
            % relative deviation of orientation and position of camera views
            deviation = max(max(sqrt(sum((Qp-Q(:,ones(1,n))).^2)),[],2), max(sqrt(sum((Tp-Tcc(:,pp)*ones(1,n)).^2)),[],2)/norm(Tcc(:,pp)));
            fprintf(1,'For view %d, relative extrinsic deviation = %1.0e;\n', pp, deviation);
            if deviation>bigeps,
                fprintf(1,'The value is larger than %1.0e, please stop and check variables ''Qp'' and ''Tp''.\n', bigeps);
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

%------------------------------------------ Multiview Optimization:
fprintf(1,'\nMultiview optimization procedure - Number of images: %d\n',length(ind_active));
fprintf(1,'Sparse Levenberg-Marquardt iterations: ');

%%% Initialization of the global parameter vector:
nview7 = (n_cam+n_ima-1)*7;
intrinsic_param = [fc;cc;alpha_c;kc];
extrinsic_param= reshape([Qcc(:, 2:n_cam), Qcw; Tcc(:, 2:n_cam), Tcw],nview7,1);

intr_update = intrinsic_param;
extr_update = extrinsic_param;
ex = []; % Global error vector
for kk =  ind_active_views,
    Qw1 = Qcw(:, kk);              % Rc1 = trans_quat_mat(Qw1);
    Tw1 = Tcw(:, kk);
    jj = (kk-1)*n_cam+1;
    X_kk = X_cell{jj};
    x_kk = x_cell{jj};
    x = project_points_mirror(X_kk,Qw1,Tw1,hand1,fc,cc,kc,alpha_c);
    ex_kk = x_kk - x;
    ex = [ex, ex_kk];
    active_view = find(active_imgviews(:, kk)');
    for pp = active_view,
        % Reproject the patterns on the images, and compute the pixel errors:
        Qp = Qcc(:,pp);
        Tp = Tcc(:,pp);
        Qw1t = Qw1;
        handkk = handcc(pp);
        if handkk ~= 1,
            Qw1t(1:2) = -Qw1t(1:2);      % Qw1t = trans_quat_mat(Hzkk*Rc1*Hzkk)
        end;
        Qwkk = quatmul(Qp, Qw1t);
        Twkk = rigid_trans(Tw1,Qp,Tp,handkk);

        kth = jj+pp-1;
        X_kk = X_cell{kth};
        x_kk = x_cell{kth};
        handkk = handkk*hand1;
        x = project_points_mirror(X_kk,Qwkk,Twkk,handkk,fc,cc,kc,alpha_c);
        ex_kk = x_kk - x;
        ex = [ex, ex_kk];
    end;
end;

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
        U = sparse([],[],[],10,10,100);
        V = sparse([],[],[],nview7,nview7,7*nview7);
        W = sparse([],[],[],10,nview7,10*nview7);

        f = intrinsic_param(1:2);
        c = intrinsic_param(3:4);
        alpha = intrinsic_param(5);
        k = intrinsic_param(6:10);

        ea = zeros(10,1);               % A'*ex
        eb = zeros(nview7,1);         % B'*ex
        for kk = ind_active_views,
            ii = 7*(n_cam+kk-2);
            Qw1 = extrinsic_param(ii+1 : ii+4);                     % Rc1 = trans_quat_mat(Qw1);
            Tw1 = extrinsic_param(ii+5 : ii+7);
            kth = (kk-1)*n_cam+1;
            X_kk = X_cell{kth};
            x_kk = x_cell{kth};
            if est_aspect_ratio,
                [x,dxdQ1,dxdT1,dxdf,dxdc,dxdk,dxdalpha] = project_points_mirror(X_kk,Qw1,Tw1,hand1,f,c,k,alpha);
            else
                [x,dxdQ1,dxdT1,dxdf,dxdc,dxdk,dxdalpha] = project_points_mirror(X_kk,Qw1,Tw1,hand1,f(1),c,k,alpha);
                dxdf = repmat(dxdf,[1 2]);
            end;
            ex_kk = x_kk - x;
            Akk = [dxdf dxdc dxdalpha dxdk];
            Bkk = [dxdQ1 dxdT1];

            U = U + sparse(Akk'*Akk);
            ea = ea + Akk'*ex_kk(:);
            % Check if this view is ill-conditioned:
            if check_cond && (cond(Bkk)> thresh_cond),
                active_imgviews(1, kk) = 0;
                fprintf(1,'\nWarning: View 1 of frame %d ill-conditioned. This view is now set inactive. \n\n',kk);
                extrinsic_param(ii+1 : ii+7) = NaN(7,1);
            else
                V(ii+1 : ii+7, ii+1 : ii+7) = V(ii+1 : ii+7, ii+1 : ii+7) + sparse(Bkk'*Bkk);
                W(1:10, ii+1 : ii+7) = W(1:10, ii+1 : ii+7) + sparse(Akk'*Bkk);
                eb(ii+1 : ii+7) = eb(ii+1 : ii+7) + Bkk'*ex_kk(:);
                active_view = find(active_imgviews(:, kk)');
                for pp = active_view,
                    % Reproject the patterns on the images, and compute the pixel errors:
                    jj = 7*(pp-2);
                    Qp = extrinsic_param(jj+1 : jj+4);
                    Tp = extrinsic_param(jj+5 : jj+7);
                    Qw1t = Qw1;
                    handkk = handcc(pp);
                    if handkk ~= 1,
                        Qw1t(1:2) = -Qw1t(1:2);      % Qw1t = trans_quat_mat(Hzkk*Rc1*Hzkk)
                    end;
                    [Qwkk,dQwkkdQp,dQwkkdQ1] = quatmul(Qp, Qw1t);
                    [Twkk,dTwkkdQp,dTwkkdTp,dTwkkdT1] = rigid_trans(Tw1,Qp,Tp,handkk);
                    if handkk ~= 1,
                        dQwkkdQ1(:, 1:2) = -dQwkkdQ1(:, 1:2);
                    end;
                    kth = (kk-1)*n_cam+pp;
                    X_kk = X_cell{kth};
                    x_kk = x_cell{kth};
                    handkk = handkk*hand1;
                    if est_aspect_ratio,
                        [x,dxdQwkk,dxdTwkk,dxdf,dxdc,dxdk,dxdalpha] = project_points_mirror(X_kk,Qwkk,Twkk,handkk,f,c,k,alpha);
                    else
                        [x,dxdQwkk,dxdTwkk,dxdf,dxdc,dxdk,dxdalpha] = project_points_mirror(X_kk,Qwkk,Twkk,handkk,f(1),c,k,alpha);
                        dxdf = repmat(dxdf,[1 2]);
                    end;
                    ex_kk = x_kk - x;

                    dxdQ1 = dxdQwkk*dQwkkdQ1;
                    dxdT1 = dxdTwkk*dTwkkdT1;
                    dxdQp = dxdQwkk*dQwkkdQp + dxdTwkk*dTwkkdQp;
                    dxdTp = dxdTwkk*dTwkkdTp;

                    Akk = [dxdf dxdc dxdalpha dxdk];
                    Bkk = [dxdQ1 dxdT1];
                    Bcc = [dxdQp dxdTp];

                    U = U + sparse(Akk'*Akk);
                    ea = ea + Akk'*ex_kk(:);
                    % Check if this view is ill-conditioned:
                    if check_cond && ((cond(Bkk)> thresh_cond) || (cond(Bcc)> thresh_cond)),
                        active_imgviews(pp, kk) = 0;
                        fprintf(1,'\nWarning: View %d of frame %d ill-conditioned. This view is now set inactive. \n\n',pp,kk);
                    else
                        V(ii+1 : ii+7, ii+1 : ii+7) = V(ii+1 : ii+7, ii+1 : ii+7) + sparse(Bkk'*Bkk);
                        V(jj+1 : jj+7, jj+1 : jj+7) = V(jj+1 : jj+7, jj+1 : jj+7) + sparse(Bcc'*Bcc);
                        W(1:10, ii+1 : ii+7) = W(1:10, ii+1 : ii+7) + sparse(Akk'*Bkk);
                        W(1:10, jj+1 : jj+7) = W(1:10, jj+1 : jj+7) + sparse(Akk'*Bcc);
                        eb(ii+1 : ii+7) = eb(ii+1 : ii+7) + Bkk'*ex_kk(:);
                        eb(jj+1 : jj+7) = eb(jj+1 : jj+7) + Bcc'*ex_kk(:);
                    end;
                end;
            end;
        end;
        active_view1 = active_imgviews(1, :);
        active_imgviews = active_imgviews.*active_view1(ones(n_cam,1),:);
        ind_active_views = find(active_view1);
        active_view = any(active_imgviews,2)';
        active_view = [active_view(2:n_cam), active_view1];

        % The following vector helps to select the variables to update (for only active images):
        selected_invar = [est_fc;center_optim*ones(2,1);est_alpha;est_dist]';
        selected_exvar = reshape(ones(7,1)*active_view,1,nview7);
        if ~est_aspect_ratio && est_fc(1),
            selected_invar(2) = 0;
        end;

        ind_va = find(selected_invar);
        ind_vb = find(selected_exvar);

        U = U(ind_va,ind_va);
        V = V(ind_vb,ind_vb);
        W = W(ind_va,ind_vb);

        ea = ea(ind_va);
        eb = eb(ind_vb);
    end;

    U_lm = U + diag(lamda*diag(U));  % U + lamda*speye(size(U));
    V_lm =V + diag(lamda*diag(V));    % V + lamda*speye(size(V));
    Y = W/V_lm;

    intr_innov = (U_lm-Y*W')\(ea-Y*eb);                   % da
    extr_innov = V_lm\(eb-W'*intr_innov);                  % db
    intr_update(ind_va) = intrinsic_param(ind_va) + intr_innov;     % updated parameters
    extr_update(ind_vb) = extrinsic_param(ind_vb) + extr_innov;

    % New intrinsic parameters:
    if ~est_aspect_ratio && all(est_fc),
        intr_update(2) = intr_update(1);
    end;
    f = intr_update(1:2);
    c = intr_update(3:4);
    alpha = intr_update(5);
    k = intr_update(6:10);
    if center_optim && (c(1)<0 || c(1)>nx || c(2)<0 || c(2)>ny),
        fprintf(1,'\nWarning: it appears that the principal point cannot be estimated. Setting center_optim = 0\n');
        center_optim = 0;
        c = cc;
        intr_update(3:4) = cc;
    end;

    % compute reprojection error vector
    ex = [];
    for kk = ind_active_views,
        ii = 7*(n_cam+kk-2);
        Qw1 = extr_update(ii+1 : ii+4);
        Qw1 = Qw1/norm(Qw1);
        extr_update(ii+1 : ii+4) = Qw1;
        Tw1 = extr_update(ii+5 : ii+7);
        kth = (kk-1)*n_cam+1;
        X_kk = X_cell{kth};
        x_kk = x_cell{kth};
        x = project_points_mirror(X_kk,Qw1,Tw1,hand1,f,c,k,alpha);
        ex_kk = x_kk - x;
        ex = [ex, ex_kk];
        active_view = find(active_imgviews(:, kk)');
        for pp = active_view,
            jj = 7*(pp-2);
            Qp = extr_update(jj+1 : jj+4);
            Qp = Qp/norm(Qp);
            extr_update(jj+1 : jj+4) = Qp;
            Tp = extr_update(jj+5 : jj+7);
            Qw1t = Qw1;
            handkk = handcc(pp);
            if handkk ~= 1,
                Qw1t(1:2) = -Qw1t(1:2);      % Qw1t = trans_quat_mat(Hzkk*Rc1*Hzkk)
            end;
            Qwkk = quatmul(Qp, Qw1t);
            Twkk = rigid_trans(Tw1,Qp,Tp,handkk);
            kth = (kk-1)*n_cam+pp;
            X_kk = X_cell{kth};
            x_kk = x_cell{kth};
            handkk = handkk*hand1;
            x = project_points_mirror(X_kk,Qwkk,Twkk,handkk,f,c,k,alpha);
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
    fprintf(1,'\nWARNING: Cannot calibrate all views of image %d.\n',kk);
    fprintf(1,'Set active_images(%d)=0;\n',kk);
end;
if ~isempty(desactivated_images),
    fprintf(1,['List of images left desactivated: ' num2str(desactivated_images) '.\n']);
end;
active_images = any(active_imgviews,1);
ind_active_views = find(active_imgviews(1, :));
ind_active = find(active_images);

%%%--------------------------- Computation of the error of estimation:

fprintf(1,'\nEstimation of uncertainties...');

% Extraction of the paramters for computing the reprojection error:
fc = intrinsic_param(1:2);
cc = intrinsic_param(3:4);
alpha_c = intrinsic_param(5);
kc = intrinsic_param(6:10);

nn = n_ima+n_cam-1;
extrinsic_param = reshape(extrinsic_param, 7, nn);
Qcc(:,2:n_cam) = extrinsic_param(1:4, 1:n_cam-1);
Tcc(:,2:n_cam) = extrinsic_param(5:7, 1:n_cam-1);
Qcw = extrinsic_param(1:4, n_cam:nn);
Tcw = extrinsic_param(5:7, n_cam:nn);

% Reproject the patterns on the images, and compute the pixel errors:
y_cell = cell(1, n_view);  % Reprojected points
ex_cell = y_cell;             % Reprojected error
ex = [];
for kk = ind_active_views,
    Qw1 = Qcw(:,kk);
    Tw1 = Tcw(:,kk);
    kth = (kk-1)*n_cam+1;
    X_kk = X_cell{kth};
    x_kk = x_cell{kth};
    y_kk = project_points_mirror(X_kk,Qw1,Tw1,hand1,fc,cc,kc,alpha_c);
    ex_kk = x_kk - y_kk;
    ex = [ex, ex_kk];
    y_cell{kth} = y_kk;
    ex_cell{kth} = ex_kk;
    active_view = find(active_imgviews(:, kk)');
    for pp = active_view,
        Qp = Qcc(:,pp);
        Tp = Tcc(:,pp);
        Qw1t = Qw1;
        handkk = handcc(pp);
        if handkk ~= 1,
            Qw1t(1:2) = -Qw1t(1:2);      % Qw1t = trans_quat_mat(Hzkk*Rc1*Hzkk)
        end;
        Qwkk = quatmul(Qp, Qw1t);
        Twkk = rigid_trans(Tw1,Qp,Tp,handkk);
        kth = (kk-1)*n_cam+pp;
        X_kk = X_cell{kth};
        x_kk = x_cell{kth};
        handkk = handkk*hand1;
        y_kk = project_points_mirror(X_kk,Qwkk,Twkk,handkk,fc,cc,kc,alpha_c);
        ex_kk = x_kk - y_kk;
        ex = [ex, ex_kk];
        y_cell{kth} = y_kk;
        ex_cell{kth} = ex_kk;
    end;
end;
err_std = std(ex,0,2);
sigma_x = std(ex(:));

% Compute the the standard deviation of parameters from std(ex):
% ex = X-f(P),  cov(param,param) = inv((JJ'* inv(cov(ex,ex))* JJ))
JJ2_inv = inv([U, W; W', V]); % not bad for sparse matrices!!
param_error = zeros(nview7+10,1);          % take value as perfect if not optimized (no error).
param_error([ind_va,ind_vb+10]) =  3*sqrt(abs(full(diag(JJ2_inv))))*sigma_x;     % 3 sigma principle
if ~est_aspect_ratio && all(est_fc),
    param_error(2) = param_error(1);
end;

%%% Extraction of the final intrinsic and extrinsic paramaters:
fc_error = param_error(1:2);
cc_error = param_error(3:4);
alpha_c_error = param_error(5);
kc_error = param_error(6:10);

param_error = reshape(param_error(11 : nview7+10),7,nn);
Qcc_error(:,2:n_cam) = param_error(1:4, 1:n_cam-1);
Tcc_error(:,2:n_cam) = param_error(5:7, 1:n_cam-1);
Qcw_error = param_error(1:4, n_cam:nn);
Tcw_error = param_error(5:7, n_cam:nn);

fprintf(1,'Done with the uncertainties!\n');

err_cam = zeros(2,n_cam);
fprintf(1,'\n\nFinal calibration results (with uncertainties):\n\n');
fprintf(1,'Focal Length:        fc = [ %3.5f   %3.5f ] ? [ %3.5f   %3.5f ]\n',[fc;fc_error]);
fprintf(1,'Principal point:     cc = [ %3.5f   %3.5f ] ? [ %3.5f   %3.5f ]\n',[cc;cc_error]);
fprintf(1,'Skew:           alpha_c = %3.5f  ?  %3.5f  =>  angle of pixel axes = %3.5f  ?  %3.5f degrees\n', [alpha_c; alpha_c_error; 90-atan(alpha_c)*180/pi; atan(alpha_c_error)*180/pi]);
fprintf(1,'Distortion:            kc = [ %3.5f   %3.5f   %3.5f   %3.5f   %3.5f ] ? [ %3.5f   %3.5f   %3.5f   %3.5f   %3.5f ]\n',[kc;kc_error]);
fprintf(1,'Pixel error of every camera view:\n');
for pp = 1:n_cam,
    err_cam(:,pp) = std(cell2mat(ex_cell(pp:n_cam:end)),0,2);
    fprintf(1,'   view %d:         err = [%3.5f, %3.5f]\n',pp,err_cam(:,pp));
end;
fprintf(1,'\nNote: The pixel error are approximately three times the standard deviations (for reference).\n\n\n');
