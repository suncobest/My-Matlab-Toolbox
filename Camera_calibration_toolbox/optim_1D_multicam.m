% intrinsic and extrinsic parameters
intr_update = reshape([fc_mat; cc_mat; alpha_vec; kc_mat; Omcc; Tcc],ncam16,1);
extr_update = reshape([Xo; thph],nima5,1);
if ~exist('init_param','var') || length(init_param) ~= ncam16+nima5,
    init_param = [intr_update; extr_update];
end;
intr_param = intr_update;
extr_param = extr_update;

% The following vector helps to select the variables to update:
selected_invar = [est_fc_mat; ones(2,1)*center_optim_vec; est_alpha_vec; est_dist_mat; active_cam(ones(6,1),:)];
selected_invar(2,:) = selected_invar(2,:).*(est_aspect_ratio_vec | ~est_fc_mat(1,:));
selected_invar(11:16,idm) = 0;
selected_exvar = reshape(active_images(ones(5,1),:),1,nima5);
ind_va = find(selected_invar);
ind_vb = find(selected_exvar);

% initial error before bundle adjustment
Xp = reshape(extr_param,5,n_ima);
Xp = gen_1D_points(Xp(1:3,:),Xp(4:5,:),rodlen);
ex = []; % Global error vector
for pp = ind_cam,
    % load camera parameters
    fc = fc_mat(:,pp);
    cc = cc_mat(:,pp);
    alpha_c = alpha_vec(pp);
    kc = kc_mat(:,pp);
    omwkk = Omcc(:,pp);
    Twkk = Tcc(:,pp);
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
            Bkk = sparse([],[],[],nb,5*na);
            idx = reshape(repmat(ind,[3,1]),1,npts3);
            id3 = reshape(repmat(active_view,[3,1]),1,nima3);
            id2 = reshape(repmat(active_view,[2,1]),1,nima2);
            id5 = reshape(repmat(active_view,[5,1]),1,nima5);
            dxdXo = dxdXp*dXpdXo(idx,id3);
            dxdtp = dxdXp*dXpdtp(idx,id2);
            for jj = 1:na,
                Bkk(:,(jj-1)*5+1:jj*5) = [dxdXo(:,(jj-1)*3+1:jj*3),dxdtp(:,(jj-1)*2+1:jj*2)];
            end;
            % Bkk = reshape([reshape(dxdXp*dXpdXo(idx,id3),[nb,3,na]), reshape(dxdXp*dXpdtp(idx,id2),[nb,2,na])], nb,5*na);

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

%%%--------------------------- Computation of the error of estimation:

% fprintf(1,'\nEstimation of uncertainties...\n');

% Extraction of the paramters for computing the reprojection error:
intr_param = reshape(intr_param, 16, n_cam);
fc_mat = intr_param(1:2,:);
cc_mat = intr_param(3:4,:);
alpha_vec = intr_param(5,:);
kc_mat = intr_param(6:10,:);
Omcc = intr_param(11:13,:);
Tcc = intr_param(14:16,:);

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
    omwkk = Omcc(:,pp);
    Twkk = Tcc(:,pp);
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
[~,ind] = max(sum(ex.^2,1));
ex_max = ex(:,ind);

return;

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
