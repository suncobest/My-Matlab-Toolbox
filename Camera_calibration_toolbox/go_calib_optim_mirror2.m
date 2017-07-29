%go_calib_optim_mirror2
%
% version with pure LM algorithm
% see go_calib_optim_mirror using quaternion rotation
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
%        Omw_mat: axis angle rotation vectors attached to the grid positions in space
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
    check_cond = 1; % Set this variable to 0 in case you don't want to extract view dynamically
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
    FOV_angle = 35; %field of view in degrees: for135 camera, 35 degree of FOV is about 57 mm focal length
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
Omw_mat = NaN(3, n_view);
Tw_mat = Omw_mat;
for pp = 1:n_cam,
    handkk = hand_list(pp);
    active_view = find(active_imgviews(pp,:));
    for kk = active_view,
        kth = (kk-1)*n_cam+pp;
        x_kk = x_cell{kth};
        X_kk = X_cell{kth};
        if isempty(x_kk) || isnan(x_kk(1)),   % check x_cell
            fprintf(1,'Warning: Cannot calibrate with view %d of image %d. This view is now set inactive.\n',pp,kk);
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
                fprintf(1,'\nWarning: View %d of image %d ill-conditioned. This view is now set inactive.\n',pp,kk);
            elseif any(isnan(omwkk)),
                active_imgviews(pp,kk) = 0;
                fprintf(1,'\nWarning: The extrinsic rotation is NaN for view %d of image %d. The view is now set inactive.\n',pp,kk);
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
nview6 = 6*n_view;
intrinsic_param = [fc;cc;alpha_c;kc];
extrinsic_param = reshape([Omw_mat; Tw_mat],nview6,1);

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
        f = intrinsic_param(1:2);
        c = intrinsic_param(3:4);
        alpha = intrinsic_param(5);
        k = intrinsic_param(6:10);

        % JJ2 = JJ'*JJ = [U, W; W', V]
        U = sparse([],[],[],10,10,100);
        V = sparse([],[],[],nview6,nview6,6*nview6);
        W = sparse([],[],[],10,nview6,10*nview6);
        ea = zeros(10,1);            % A'*ex
        eb = zeros(nview6,1);      % B'*ex

        for kk = ind_active_views,
            ii = 6*(kk-1);
            omwkk = extrinsic_param(ii+1 : ii+3);
            Twkk = extrinsic_param(ii+4 : ii+6);
            if any(isnan(omwkk)),
                fprintf(1,'Extrinsic parameters at view %d of frame %d do not exist!\n',mod(kk-1,n_cam)+1,ceil(kk/n_cam));
                return;
            end;
            X_kk = X_cell{kk};
            x_kk = x_cell{kk};
            handkk = hand_list(mod(kk-1,n_cam)+1);
            if est_aspect_ratio,
                [x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points_mirror2(X_kk,omwkk,Twkk,handkk,f,c,k,alpha);
            else
                [x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points_mirror2(X_kk,omwkk,Twkk,handkk,f(1),c,k,alpha);
                dxdf = repmat(dxdf,[1 2]);
            end;
            ex_kk = x_kk - x;
            Akk = [dxdf dxdc dxdalpha dxdk];
            Bkk = [dxdom dxdT];

            U = U + sparse(Akk'*Akk);
            ea = ea + Akk'*ex_kk(:);
            % Check if this view is ill-conditioned:
            if check_cond && (cond(Bkk)> thresh_cond),
                fprintf(1,['\nWarning: View %d of frame %d ill-conditioned. This view is now set inactive. \n' ...
                    '(note: to disactivate this option, set check_cond=0)\n'],mod(kk-1,n_cam)+1,ceil(kk/n_cam));
                active_imgviews(kk) = 0;
                extrinsic_param(ii+1 : ii+6) = NaN(6,1);
            else
                V(ii+1 : ii+6, ii+1 : ii+6) = sparse(Bkk'*Bkk);
                W(1:10, ii+1 : ii+6) = sparse(Akk'*Bkk);
                eb(ii+1 : ii+6) = Bkk'*ex_kk(:);
            end;
        end;
        active_view = active_imgviews(:)';
        ind_active_views = find(active_view);

        % The following vector helps to select the variables to update (for only active images):
        selected_invar = [est_fc; center_optim*ones(2,1); est_alpha; est_dist]';
        selected_exvar = reshape(ones(6,1)*active_view,1,nview6);
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
        ii = 6*(kk-1);
        omwkk = extr_update(ii+1 : ii+3);
        Twkk = extr_update(ii+4 : ii+6);
        X_kk = X_cell{kk};
        x_kk = x_cell{kk};
        handkk = hand_list(mod(kk-1,n_cam)+1);
        x = project_points_mirror2(X_kk,omwkk,Twkk,handkk,f,c,k,alpha);
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
    fprintf(1,['List of images left desactivated: ' num2str(desactivated_images) '.\n' ] );
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
extrinsic_param = reshape(extrinsic_param, 6, n_view);
Omw_mat = extrinsic_param(1:3,:);
Tw_mat = extrinsic_param(4:6,:);

y_cell = cell(1, n_view);  % Reprojected points
ex_cell = y_cell;             % Reprojected error
H_cell = y_cell;            % recomputer the collineations
ex = [];
ind_active_views = find(active_imgviews(:)');
for kk = ind_active_views,
    omwkk = Omw_mat(:, kk);
    Twkk = Tw_mat(:, kk);
    X_kk = X_cell{kk};
    x_kk = x_cell{kk};
    handkk = hand_list(mod(kk-1,n_cam)+1);
    y_kk = project_points_mirror2(X_kk,omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
    ex_kk = x_kk-y_kk;
    ex = [ex, ex_kk];
    y_cell{kk} = y_kk;
    ex_cell{kk} = ex_kk;
    Rwkk = rodrigues(omwkk);
    Hkk = KK * [Rwkk(:,1) Rwkk(:,2) Twkk];
    H_cell{kk} = Hkk / Hkk(3,3);
end;

err_std = std(ex,0,2);
sigma_x = std(ex(:));

% Compute the the standard deviation of parameters from std(ex):
% ex = X-f(P),  cov(param,param) = inv((JJ'* inv(cov(ex,ex))* JJ))
JJ2 = [U, W; W', V]; % not bad for sparse matrices!!
JJ2_inv = speye(size(JJ2))/JJ2;
param_error = zeros(nview6+10,1);        % take value as perfect if not optimized (no error).
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
param_error = reshape(param_error(11 : nview6+10),6,n_view);
Omw_mat_error = param_error(1:3,:);
Tw_mat_error = param_error(4:6,:);

fprintf(1,'Done with the uncertainties!\n');

err_cam = zeros(2,n_cam);
fprintf(1,'\n\nCalibration results after optimization (with uncertainties):\n\n');
fprintf(1,'Focal Length:      fc = [%3.5f, %3.5f] ? [%3.5f, %3.5f]\n',[fc;fc_error]);
fprintf(1,'Principal point:   cc = [%3.5f, %3.5f] ? [%3.5f, %3.5f]\n',[cc;cc_error]);
fprintf(1,'Skew:         alpha_c = %3.5f ? %3.5f => Skew angle = %3.5f ? %3.5f degrees\n', ...
    [alpha_c; alpha_c_error; 90-atan(alpha_c)*180/pi; atan(alpha_c_error)*180/pi]);
fprintf(1,'Distortion:          kc = [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f] ? [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f]\n',[kc;kc_error]);
fprintf(1,'Pixel error of every camera view:\n');
for pp = 1:n_cam,
    err_cam(:,pp) = std(cell2mat(ex_cell(pp:n_cam:end)),0,2);
    fprintf(1,'   view %d:         err = [%3.5f, %3.5f]\n',pp,err_cam(:,pp));
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
            fprintf(1,'No transformation of view %d wrt view 1 found!\n',pp);
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
            % relative deviation of orientation and position of camera views
            deviation = max(max(sqrt(sum((Qck-Q(:,ones(1,n))).^2)),[],2), max(sqrt(sum((Tck-Tcc(:,pp)*ones(1,n)).^2)),[],2)/norm(Tcc(:,pp)));
            fprintf(1,'For view %d, relative extrinsic deviation = %1.0e;\n', pp, deviation);
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
            fprintf(1,'Image %d is inactive!\n',kk);
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
    optim_split_extrinsic;
end;
flag = input('Further refine intrinsic parameters or not? ([]=no, other=yes)','s');
if isempty(flag)
    disp('Done.');
    return;
end;

for count=1:3,
    fprintf(1,'\nRefine intrinsic and extrinsic parameters: %2d\n', count);
     optim_split_views;
end;
flag = input('Change some estimation settings to further refine parameters? ([]=no, other=yes)','s');
if ~isempty(flag),
    fprintf(1,'\nCarefully choose to change estimation settings, especially distortion settings!\n');
    if ~isequal(est_fc,[1;1]),
        if isequal(est_fc,[1;0]),
            fprintf(1,'The second component of focal (fc(2)) is not estimated.\n');
        else
            if isequal(est_fc,[0;1]),
                fprintf(1,'The first component of focal (fc(1)) is not estimated.\n');
            else
                fprintf(1,'The focal vector fc is not optimized.\n');
            end;
        end;
        flag = input('Do you want to change it? ([]=no, other=yes)','s');
        if ~isempty(flag),
            est_fc = [1;1];
        end;
    end;

    if ~est_aspect_ratio,
        fprintf(1,'The aspect ratio of pixels is not estimated (fc(1)==fc(2)).\n');
        flag = input('Do you want to estimate it? ([]=no, other=yes)','s');
        if ~isempty(flag),
            est_aspect_ratio = 1;
        end;
    end;

    if ~center_optim,
        fprintf(1,'The principal point cc is not optimized.\n');
        flag = input('Do you want to estimate it? ([]=no, other=yes)','s');
        if ~isempty(flag),
            center_optim = 1;
        end;
    end;

     if ~est_alpha,
        fprintf(1,'The skew of pixels is not optimized (alpha_c=0).\n');
        flag = input('Do you want to estimate it? ([]=no, other=yes)','s');
        if ~isempty(flag),
            est_alpha = 1;
        end;
    end;

    if ~norm(est_dist),
        fprintf(1,'\nThe distortion coefficients are not estimated.\n');
        flag = input('Do you want to estimate it? ([]=no, other=yes with cautiousness)','s');
        if ~isempty(flag),
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
    end;

    for count=1:2,
        fprintf(1,'\nFurther refine intrinsic and extrinsic parameters: %2d\n', count);
        optim_split_views;
    end;
    fprintf(1,'\nRefine extrinsic parameters in the end:\n');
    optim_multicams_extrinsic;
end;
disp('Done.');
