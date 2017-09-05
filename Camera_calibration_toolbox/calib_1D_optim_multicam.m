% calib_1D_optim_multicam
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
            clear Xrod Xori Xdir fc_mat cc_mat kc_mat alpha_vec Omcw Tcw;
        end;
    else
        fprintf(1,'\nThere is no data required for calibration!\n');
        return;
    end;
end;

active_images = sum(active_imgviews,1)>1;
ind_active = find(active_images);
active_imgviews(:,~active_images) = 0;

if ~exist('MaxIter','var'),
    MaxIter = 30; % Maximum number of iterations in the main LM algorithm
end;

if exist('fc_mat','var') && size(fc_mat,2)==n_cam,
    fprintf(1,'\nIt seems that the intrinsic parameters exist!\nIntrinsic parameters do not need estimation if they are perfect.\n');
    est_intrinsic = input('Refine intrinsic parameters or not? ([]=no, other=yes) ','s');
    est_intrinsic = ~isempty(est_intrinsic);
else
    est_intrinsic = 1;
end;

if est_intrinsic,
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
        est_alpha_vec = input(['est_alpha_vec = ([] = [' num2str(zeros(1,n_cam)) ']) ']);
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
        est_aspect_ratio_vec = input(['est_aspect_ratio_vec = ([] = [' num2str(ones(1,n_cam)) ']) ']);
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
        center_optim_vec = input(['center_optim_vec = ([] = [' num2str(ones(1,n_cam)) ']) ']);
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
    while flag,
        fprintf(1,'\nDo you want to estimate the focal length of all %d cameras?\nSet to zero if you don''t!\n',n_cam);
        % by default estimate focal length
        est_fc_mat = input(['est_fc_mat = [1;1] * ([] = [' num2str(ones(1,n_cam)) ']) ']);
        if isempty(est_fc_mat),
            est_fc_mat = ones(2,n_cam);
            flag = 0;
        else
            est_fc_mat = est_fc_mat(:)';
            flag = length(est_fc_mat)~=n_cam;
            if flag,
                fprintf(1,'\nDimension of est_fc_mat do not match number of cameras! Please input again!\n');
            else
                est_fc_mat = est_fc_mat(ones(2,1),:);
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
    while flag,
        flag = input('Number of distortion coefficients to estimate? (0:5, []=0) ');
        if isempty(flag),
            est_dist_mat = zeros(5,n_cam);
            flag = 0;
        else
            flag = round(flag);
            if flag<0 || flag>5,
                fprintf(1,'\nThe number of distortion coefficients for estimation is assumed to be 0:5!\n');
                flag = 1;
            else
                while 1,
                    est_dist = zeros(5,1);
                    est_dist(1:flag) = 1;
                    flag = sprintf('[%d;%d;%d;%d;%d]',est_dist);
                    fprintf(1,'\nDo you want to estimate the distortion of all %d cameras?\nSet to zero if you don''t!\n',n_cam);
                    est_dist_mat = input(['est_dist_mat = ' flag ' * ([] = [' num2str(ones(1,n_cam)) ']) ']);
                    if isempty(est_dist_mat),
                        est_dist_mat = est_dist(:,ones(1,n_cam));
                        break;
                    else
                        est_dist_mat = est_dist_mat(:)';
                        if length(est_dist_mat)~=n_cam,
                            fprintf(1,'\nDimension of est_dist_mat do not match number of cameras! Please input again!\n');
                        else
                            est_dist_mat = est_dist .* est_dist_mat;
                            break;
                        end;
                    end;
                end;
                flag = 0;
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
        FOV_angle = 70; %field of view in degrees: for 135 camera, 70 degree of FOV is about 25 mm focal length。
        fprintf(1,'Initialization of the focal length with FOV of %3.1f degrees.\n\n',FOV_angle);
        fc_mat = ones(2,1)*(imsize(1,:)/2)/tan(pi*FOV_angle/360);    % FOV_angle=2*atan(nx/(2*fc))

        % flag = input('The calibration rod was only under rotation or not? ([]=no, other=yes) ','s');
        % if isempty(flag),
        %     FOV_angle = 70; %field of view in degrees: for 135 camera, 70 degree of FOV is about 25 mm focal length。
        %     fprintf(1,'Initialization of the focal length with FOV of %3.1f degrees.\n\n',FOV_angle);
        %     fc_mat = ones(2,1)*(imsize(1,:)/2)/tan(pi*FOV_angle/360);    % FOV_angle=2*atan(nx/(2*fc))
        % else
        %     fprintf(1,'\nInitialization of the intrinsic parameters using Zhang Zhengyou''s algorithm.\n');
        %     intrinsic_1D_rotation;
        % end;
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
        fprintf(1,'\nInitialization of all camera image skew to zero.\n');
        alpha_vec = zeros(1,n_cam);
    end;

    % no principal point estimation, no skew estimation
    est_alpha_vec(center_optim_vec==0) = 0;

    % alpha = 0 when skew is not estimated
    alpha_vec(est_alpha_vec==0) = 0;

    % set to zero the distortion coefficients that are not estimated
    kc_mat = kc_mat .* est_dist_mat;
end;

% threshold to terminate the main LM iteration
gradeps = 1e-5;

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
                if est_intrinsic,
                    % refine camera pair
                    est_fc = est_fc_mat(:,id);
                    center_optim = center_optim_vec(id);
                    est_dist = est_dist_mat(:,id);
                    est_alpha = est_alpha_vec(id);
                    est_aspect = est_aspect_ratio_vec(id);
                    [~,om2,T2,fc,cc,kc,alpha_c] = binocular_1D_optim(xx,rodlen,om2,T2,handkk,fc,cc,kc,alpha_c,...
                                                                      est_fc,center_optim,est_dist,est_alpha,est_aspect);
                    [om,T] = compose_motion2(Omcc(:,id(1)),Tcc(:,id(1)),om2,T2,handkk);
                    fc_mat(:,id) = fc;
                    cc_mat(:,id) = cc;
                    kc_mat(:,id) = kc;
                    alpha_vec(id) = alpha_c;
                else
                    % triangulation and  determine the scale factor
                    XX = compute_structure2(xx,[zeros(3,1),om2],[zeros(3,1),T2],[1,handkk],fc,cc,kc,alpha_c);
                    XX = reshape(XX,[3,np1D,length(kk)]);
                    ind = all(all(~isnan(XX),1),2);
                    Xlen = permute(sqrt(sum(diff(XX(:,:,ind),[],2).^2,1)),[3,2,1]);  % length of rod (with ||T2||=1)
                    s = mean(diff(rodlen)./mean(Xlen,1));
                    [om,T] = compose_motion2(Omcc(:,id(1)),Tcc(:,id(1)),om2,T2*s,handkk);
                end;
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
    err0 = zeros(6,nc);
    for kk = 1:nc,
        pp = ind_cam(kk);
        [om,T] = compose_motion2(Omcw(:,idm),Tcw(:,idm),Omcc(:,pp),Tcc(:,pp),handcc(pp));
        Om2(:,pp) = om;
        T2(:,pp) = T;
        err0(1:3,kk) = om-Omcw(:,pp);
        err0(4:6,kk) = (T-Tcw(:,pp))/norm(Tcw(:,pp));
    end;
    if exist('Xrod','var'),
        ind = reshape(repmat(active_images,[np1D,1]),1,npts);
        errX0 = XX(:,ind)-rigid_refmotion(Xrod(:,ind),Omcw(:,idm),Tcw(:,idm),hand_list(idm));
        estdX0 = std(errX0,0,2);
    end;
end;


%% ------------------------------------------ Main Optimization:

fprintf(1,'\nMain calibration optimization procedure - Bundle adjustment with %d cameras\n', n_cam);
fprintf(1,'Sparse Levenberg-Marquardt iterations: ');

%%% Initialization of the global parameter vector:
ncam16 = 16*n_cam;
ncam6 = 6*n_cam;
nima5 = 5*n_ima;
nima3 = 3*n_ima;
nima2 = 2*n_ima;
npts3 = 3*npts;

tstart = tic;

if est_intrinsic,
    fprintf(1,'\nRefine intrinsic and extrinsic parameters:\n');
    optim_1D_multicam;
else
    fprintf(1,'\nRefine extrinsic parameters:\n');
    optim_1D_extrinsic;
end;

telapsed = toc(tstart);


% flag = input('\nFurther refine the extrinsic parameters or not? ([]=no, other=yes) ','s');
% if ~isempty(flag)
%     fprintf(1,'\nRefine extrinsic parameters in the end:\n');
%     optim_1D_extrinsic;
% end;

XX = compute_structure2(xx,Omcc,Tcc,handcc,fc_mat,cc_mat,kc_mat,alpha_vec);
ind = reshape(repmat(active_images,[np1D,1]),1,npts);
errX = Xp(:,ind)-XX(:,ind);
estdX = std(errX,0,2);
[~,ind] = max(sum(errX.^2,1));
errX_max = errX(:,ind);

save_name = 'Calib_Results_1D';
fprintf(1,['\nSaving calibration results under ' save_name '.mat\n']);
string_save = ['save ' save_name ' n_cam n_ima imsize fc_mat cc_mat kc_mat alpha_vec handcc hand_list np1D' ...
                       ' rodlen Omcc Tcc idm A_cam costs paths pathm ind_active active_images active_imgviews' ...
                       ' est_fc_mat center_optim_vec est_alpha_vec est_dist_mat est_aspect_ratio_vec Xp Xo thph' ...
                       '  x_cell y_cam ex_cam err_cam err_cam ex err_std err_std0 ex_max errX estdX errX_max'];

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
    string_save = ['save ' save_name ' Omcw Tcw'];
    if exist('Xrod','var'),
        ind = reshape(repmat(active_images,[np1D,1]),1,npts);
        errX = Xp(:,ind)-rigid_refmotion(Xrod(:,ind),Omcw(:,idm),Tcw(:,idm),hand_list(idm));
        estdX = std(errX,0,2);
        string_save = ['save ' save_name ' Xrod'];
    end;
end;

eval(string_save);
disp('Done.');
