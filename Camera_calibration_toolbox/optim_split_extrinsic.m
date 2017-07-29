%% ------------------------------------------ Multiview Optimization:
if ~exist('desactivated_images','var') || ~exist('thresh_cond','var'),
    desactivated_images = [];
    % Conditioning threshold for view rejection
    thresh_cond = 1e8;
    % threshold to terminate the main LM iteration
    gradeps = eps*1e8;
    n_view = n_ima * n_cam;
end;
fprintf(1,'\nMultiview optimization procedure - Number of images: %d\n',length(ind_active));
fprintf(1,'Sparse Levenberg-Marquardt iterations: ');

%%% Initialization of the global parameter vector:
nview6 = (n_cam+n_ima)*6;
extrinsic_param= reshape([Omcc, Omcw; Tcc, Tcw],nview6,1);

extr_update = extrinsic_param;
ex = []; % Global error vector
for pp = 1:n_cam,
    omck = Omcc(:,pp);
    Tck = Tcc(:,pp);
    hand = handcc(pp);
    handkk = hand_list(pp);
    active_view = find(active_imgviews(pp,:));
    for kk = active_view,
        % Reproject the patterns on the images, and compute the pixel errors:
        omw1 = Omcw(:, kk);              % Rc1 = rodrigues(omw1);
        Tw1 = Tcw(:, kk);
        [omwkk,Twkk] = compose_motion2(omw1,Tw1,omck,Tck,hand);
        kth = (kk-1)*n_cam+pp;
        X_kk = X_cell{kth};
        x_kk = x_cell{kth};
        x = project_points_mirror2(X_kk,omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
        ex_kk = x_kk - x;
        ex = [ex, ex_kk];
        end;
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
        % JJ2 = JJ'*JJ = V
        V = sparse([],[],[],nview6,nview6,6*nview6);
        eb = zeros(nview6,1);         % B'*ex
        for pp = 1:n_cam,
            jj = 6*(pp-1);
            omck = extrinsic_param(jj+1 : jj+3);
            Tck = extrinsic_param(jj+4 : jj+6);
            hand = handcc(pp);
            handkk = hand_list(pp);
            active_view = find(active_imgviews(pp,:));
            for kk = active_view,
                % Reproject the patterns on the images, and compute the pixel errors:
                ii = 6*(n_cam+kk-1);
                omw1 = extrinsic_param(ii+1 : ii+3);                     % Rc1 = rodrigues(omw1);
                Tw1 = extrinsic_param(ii+4 : ii+6);
                [omwkk,Twkk,domwkkdom1,domwkkdomck,dTwkkdomck,dTwkkdT1,dTwkkdTck] ...
                    = compose_motion2(omw1,Tw1,omck,Tck,hand);

                kth = (kk-1)*n_cam+pp;
                X_kk = X_cell{kth};
                x_kk = x_cell{kth};
                if est_aspect_ratio,
                    [x,dxdomwkk,dxdTwkk] = project_points_mirror2(X_kk,omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
                else
                    [x,dxdomwkk,dxdTwkk] = project_points_mirror2(X_kk,omwkk,Twkk,handkk,fc(1),cc,kc,alpha_c);
                end;
                ex_kk = x_kk - x;

                dxdom1 = dxdomwkk*domwkkdom1;
                dxdT1 = dxdTwkk*dTwkkdT1;
                dxdomck = dxdomwkk*domwkkdomck + dxdTwkk*dTwkkdomck;
                dxdTck = dxdTwkk*dTwkkdTck;
                Bkk = [dxdom1 dxdT1];
                Bcc = [dxdomck dxdTck];
                % Check if this view is ill-conditioned:
                if check_cond && ((cond(Bkk)> thresh_cond) || (cond(Bcc)> thresh_cond)),
                    active_imgviews(pp, kk) = 0;
                    fprintf(1,'\nWarning: View #%d of frame #%d ill-conditioned. This view is now set inactive. \n\n',pp,kk);
                else
                    V(ii+1 : ii+6, ii+1 : ii+6) = V(ii+1 : ii+6, ii+1 : ii+6) + sparse(Bkk'*Bkk);
                    V(jj+1 : jj+6, jj+1 : jj+6) = V(jj+1 : jj+6, jj+1 : jj+6) + sparse(Bcc'*Bcc);
                    eb(ii+1 : ii+6) = eb(ii+1 : ii+6) + Bkk'*ex_kk(:);
                    eb(jj+1 : jj+6) = eb(jj+1 : jj+6) + Bcc'*ex_kk(:);
                end;
            end;
        end;
        active_images = any(active_imgviews,1);
        active_view = [false, any(active_imgviews(2:end,:),2)', active_images];

        % The following vector helps to select the variables to update (for only active images):
        selected_exvar = reshape(ones(6,1)*active_view,1,nview6);
        ind_vb = find(selected_exvar);
        V = V(ind_vb,ind_vb);
        eb = eb(ind_vb);
    end;
    V_lm =V + diag(lamda*diag(V));    % V + lamda*speye(size(V));
    extr_innov = V_lm\eb;                  % db
    extr_update(ind_vb) = extrinsic_param(ind_vb) + extr_innov;

    % compute reprojection error vector
    ex = [];
    for pp = 1:n_cam,
        jj = 6*(pp-1);
        omck = extr_update(jj+1 : jj+3);
        Tck = extr_update(jj+4 : jj+6);
        hand = handcc(pp);
        handkk = hand_list(pp);
        active_view = find(active_imgviews(pp,:));
        for kk = active_view,
            ii = 6*(n_cam+kk-1);
            omw1 = extr_update(ii+1 : ii+3);
            Tw1 = extr_update(ii+4 : ii+6);
            [omwkk,Twkk] = compose_motion2(omw1,Tw1,omck,Tck,hand);
            kth = (kk-1)*n_cam+pp;
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
        extrinsic_param = extr_update;
        ex2 = ex2_update;
        updateJ=1;
    else
        lamda = lamda*10;
        updateJ=0;
    end;
    if norm(eb) < gradeps,
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
    fprintf(1,['List of images left desactivated: ' num2str(desactivated_images) '.\n' ] );
end;
active_images = any(active_imgviews,1);
ind_active = find(active_images);

%%%--------------------------- Computation of the error of estimation:

fprintf(1,'\nEstimation of uncertainties...');

% Extraction of the paramters for computing the reprojection error:
nn = n_ima+n_cam;
extrinsic_param = reshape(extrinsic_param, 6, nn);
Omcc(:,2:n_cam) = extrinsic_param(1:3, 2:n_cam);
Tcc(:,2:n_cam) = extrinsic_param(4:6, 2:n_cam);
Omcw = extrinsic_param(1:3, n_cam+1:nn);
Tcw = extrinsic_param(4:6, n_cam+1:nn);

% Reproject the patterns on the images, and compute the pixel errors:
y_cell = cell(1, n_view);  % Reprojected points
ex_cell = y_cell;             % Reprojected error
err_cam = zeros(2,n_cam);
ex = [];
for pp = 1:n_cam,
    omck = Omcc(:,pp);
    Tck = Tcc(:,pp);
    hand = handcc(pp);
    handkk = hand_list(pp);
    active_view = find(active_imgviews(pp,:));
    for kk = active_view,
        omw1 = Omcw(:,kk);
        Tw1 = Tcw(:,kk);
        [omwkk,Twkk] = compose_motion2(omw1,Tw1,omck,Tck,hand);
        kth = (kk-1)*n_cam+pp;
        X_kk = X_cell{kth};
        x_kk = x_cell{kth};
        y_kk = project_points_mirror2(X_kk,omwkk,Twkk,handkk,fc,cc,kc,alpha_c);
        ex_kk = x_kk - y_kk;
        ex = [ex, ex_kk];
        y_cell{kth} = y_kk;
        ex_cell{kth} = ex_kk;
    end;
    err_cam(:,pp) = std(cell2mat(ex_cell(pp:n_cam:end)),0,2);
end;
err_std = std(ex,0,2);
sigma_x = std(ex(:));

% Compute the the standard deviation of parameters from std(ex):
% ex = X-f(P),  cov(param,param) = inv((JJ'* inv(cov(ex,ex))* JJ))
JJ2_inv = speye(size(V))/V;
param_error = zeros(nview6,1);          % take value as perfect if not optimized (no error).
param_error(ind_vb) =  3*sqrt(abs(full(diag(JJ2_inv))))*sigma_x;     % 3 sigma principle

param_error = reshape(param_error,6,nn);
Omcc_error(:,2:n_cam) = param_error(1:3, 2:n_cam);
Tcc_error(:,2:n_cam) = param_error(4:6, 2:n_cam);
Omcw_error = param_error(1:3, n_cam+1:nn);
Tcw_error = param_error(4:6, n_cam+1:nn);

fprintf(1,'Done with the uncertainties!\n');
fprintf(1,'\n\nFinal calibration results (with uncertainties):\n\n');
fprintf(1,'Focal Length:      fc = [%3.5f, %3.5f] ? [%3.5f, %3.5f]\n',[fc;fc_error]);
fprintf(1,'Principal point:   cc = [%3.5f, %3.5f] ? [%3.5f, %3.5f]\n',[cc;cc_error]);
fprintf(1,'Skew:         alpha_c = %3.5f ? %3.5f => Skew angle = %3.5f ? %3.5f degrees\n', ...
    [alpha_c; alpha_c_error; 90-atan(alpha_c)*180/pi; atan(alpha_c_error)*180/pi]);
fprintf(1,'Distortion:          kc = [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f] ? [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f]\n',[kc;kc_error]);
fprintf(1,'Pixel error :        err = [%3.5f, %3.5f]\n',err_std);
fprintf(1,'Pixel error of every camera view:\n');
for pp = 1:n_cam,
    fprintf(1,'   view %d:         err = [%3.5f, %3.5f]\n',pp,err_cam(:,pp));
end;
fprintf(1,'\nNote: The pixel error are approximately three times the standard deviations (for reference).\n\n\n');
