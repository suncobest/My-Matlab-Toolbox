if ~exist('cc_mat','var')||~exist('fc_mat','var'),
    fprintf(1,'No calibration data available.\n');
    return;
end;

if ~exist('n_cam','var'),
    n_cam = 1;
end;

if ~exist('kc_mat','var'),
    kc_mat = zeros(5,n_cam);
else
    [m,n] = size(kc_mat);
    if n~=n_cam,
        fprintf(1,'\nWARNING: Dimension of distortion matrix do not match number of cameras!\n');
        kc_mat = zeros(5,n_cam);
    elseif m<5,
        kc_mat = [kc_mat;zeros(5-m,n_cam)];
    elseif m>5,
        kc_mat = kc_mat(1:5,:);
    end;
end;

if ~exist('alpha_vec','var'),
    alpha_vec = zeros(1,n_cam);
end;

if ~exist('hand_list','var'),
    hand_list = ones(1,n_cam);
end;

save_name = 'Multicam_Calib_Results';

if exist([save_name '.mat'],'file')==2,
    disp('WARNING: File ''Multicam_Calib_Results.mat'' already exists!');
    if exist('copyfile','builtin')==5,
        cont = -1;
        flag = 1;
        while flag,
            cont = cont + 1;
            save_name = ['Old_Multicam_Calib_Results' num2str(cont)];
            flag = (exist([save_name '.mat'],'file')==2);
        end;
        copyfile('Multicam_Calib_Results.mat',[save_name '.mat']);
        disp(['Copying the current ''Multicam_Calib_Results.mat'' file to ''' save_name '.mat''']);
        if exist('Multicam_Calib_Results.m','file')==2,
            copyfile('Multicam_Calib_Results.m',[save_name '.m']);
            disp(['Copying the current ''Multicam_Calib_Results.m'' file to ''' save_name '.m''']);
        end;
        cont_save = 1;
    else
        disp('The file ''Multicam_Calib_Results.mat'' is about to be changed.');
        cont_save = input('Do you want to continue? ([]=no,other=yes) ','s');
        cont_save = ~isempty(cont_save);
    end;
else
    cont_save = 1;
end;

if ~cont_save,
    fprintf(1,'\nCalibration results are not saved!\n');
    return;
end;

save_name = 'Multicam_Calib_Results';
fprintf(1,['\nSaving intrinsic calibration results under ''' save_name '.mat''\n']);
string_save = ['save ' save_name ' fc_mat cc_mat kc_mat alpha_vec n_cam hand_list imsize'];

fid = fopen([save_name '.m'],'wt');
fprintf(1,'Generating the matlab script file ''%s.m'' containing the intrinsic and extrinsic parameters...\n',save_name);
fprintf(fid,'%%%% Intrinsic and Extrinsic Camera Parameters:\n');
fprintf(fid,'%%\n');
fprintf(fid,'%% This script file can be directly excecuted under Matlab to recover the camera intrinsic and extrinsic parameters.\n');
fprintf(fid,'%% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.\n');
fprintf(fid,'%%            All those complementary variables are saved in the complete matlab data file Multicam_Calib_Results.mat.\n');
fprintf(fid,'%% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/doc/\n\n\n');

fprintf(fid,'%%-- Focal length matrix:\n');
fprintf(fid,'fc_mat = [\n');
for pp = 1:n_cam,
    fprintf(fid,'%5.15f, %5.15f;\n', fc_mat(:,pp));
end;
fprintf(fid,']'';\n\n');

fprintf(fid,'%%-- Principal point matrix:\n');
fprintf(fid,'cc_mat = [\n');
for pp = 1:n_cam,
    fprintf(fid,'%5.15f, %5.15f;\n', cc_mat(:,pp));
end;
fprintf(fid,']'';\n\n');

fprintf(fid,'%%-- Skew coefficient vector:\n');
fprintf(fid,'alpha_vec = [\n');
for pp = 1:n_cam,
    fprintf(fid,'%5.15f;\n', alpha_vec(pp));
end;
fprintf(fid,']'';\n\n');

fprintf(fid,'%%-- Distortion coefficients matrix:\n');
fprintf(fid,'kc_mat = [\n');
for pp = 1:n_cam,
    fprintf(fid,'%5.15f, %5.15f, %5.15f, %5.15f, %5.15f;\n', kc_mat(:,pp));
end;
fprintf(fid,']'';\n\n');

fprintf(fid,'%%-- Number of cameras:\n');
fprintf(fid,'n_cam = %d;\n\n',n_cam);

fprintf(fid,'%%-- Axis handness of all cameras:\n');
fprintf(fid,'hand_list = [\n');
for pp = 1:n_cam,
    fprintf(fid,'%d;\n', hand_list(pp));
end;
fprintf(fid,']'';\n\n');

fprintf(fid,'%%-- Image size matrix:\n');
fprintf(fid,'imsize = [\n');
for pp = 1:n_cam,
    fprintf(fid,'%d, %d;\n', imsize(:,pp));
end;
fprintf(fid,']'';\n\n\n');

if exist('fc_mat_error','var'),
    fprintf(1,['\nSaving intrinsic error under ' save_name '.mat\n']);
    string_save = [string_save ' fc_mat_error cc_mat_error alpha_vec_error kc_mat_error'];

    fprintf(fid,'%%-- Focal length uncertainty matrix:\n');
    fprintf(fid,'fc_mat_error = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%5.15f, %5.15f;\n', fc_mat_error(:,pp));
    end;
    fprintf(fid,']'';\n\n');

    fprintf(fid,'%%-- Principal point uncertainty matrix:\n');
    fprintf(fid,'cc_mat_error = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%5.15f, %5.15f;\n', cc_mat_error(:,pp));
    end;
    fprintf(fid,']'';\n\n');

    fprintf(fid,'%%-- Skew coefficient uncertainty vector:\n');
        fprintf(fid,'alpha_vec_error = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%5.15f;\n', alpha_vec_error(:,pp));
    end;
    fprintf(fid,']'';\n\n');

    fprintf(fid,'%%-- Distortion coefficients uncertainty matrix:\n');
    fprintf(fid,'kc_mat_error = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%5.15f, %5.15f, %5.15f, %5.15f, %5.15f;\n', kc_mat_error(:,pp));
    end;
    fprintf(fid,']'';\n\n\n');
end;

if exist('Qw_mat_error', 'var'),
    rotflag = 1;
else
    rotflag = 0;
end;

if exist('Tcc','var'),
    if ~exist('handcc','var'),
        handcc = hand_list(1)*hand_list;
    end;
    fprintf(1,'\nConfiguration of multiple cameras detected!\n');

    fprintf(fid,'%%-- Axis handness of all cameras wrt camera 1:\n');
    fprintf(fid,'handcc = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%d;\n', handcc(pp));
    end;
    fprintf(fid,']'';\n\n');

    fprintf(fid,'%%-- Refined extrinsic parameters of all cameras wrt camera 1:\n');
    if rotflag,
        string_save = [string_save ' handcc Qcc Tcc'];
        fprintf(fid,'%%-- Quaternion orientation (Qcc) and position (Tcc) vectors\n\n');
        fprintf(fid,'Qcc = [\n');
        for pp = 1:n_cam,
            fprintf(fid,'%5.15f, %5.15f, %5.15f, %5.15f;\n', Qcc(:,pp));
        end;
        fprintf(fid,']'';\n\n');
    else
        string_save = [string_save ' handcc Omcc Tcc'];
        fprintf(fid,'%%-- Axis angle orientation (Omcc) and position (Tcc) vectors\n\n');
        fprintf(fid,'Omcc = [\n');
        for pp = 1:n_cam,
            fprintf(fid,'%5.15f, %5.15f, %5.15f;\n', Omcc(:,pp));
        end;
        fprintf(fid,']'';\n\n');
    end;

    fprintf(fid,'Tcc = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%5.15f, %5.15f, %5.15f;\n', Tcc(:,pp));
    end;
    fprintf(fid,']'';\n\n');
end;

if exist('Tcc_error','var'),
    if rotflag,
        string_save = [string_save ' Qcc_error Tcc_error'];
        fprintf(fid,'%%-- Uncertainties of orientation (Qcc) and position (Tcc):\n\n');
        fprintf(fid,'Qcc_error = [\n');
        for pp = 1:n_cam,
            fprintf(fid,'%5.15f, %5.15f, %5.15f, %5.15f;\n', Qcc_error(:,pp));
        end;
        fprintf(fid,']'';\n\n');
    else
        string_save = [string_save ' Omcc_error Tcc_error'];
        fprintf(fid,'%%-- Uncertainties of orientation (Omcc) and position (Tcc):\n\n');
        fprintf(fid,'Omcc_error = [\n');
        for pp = 1:n_cam,
            fprintf(fid,'%5.15f, %5.15f, %5.15f;\n', Omcc_error(:,pp));
        end;
        fprintf(fid,']'';\n\n');
    end;

    fprintf(fid,'Tcc_error = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%5.15f, %5.15f, %5.15f;\n', Tcc_error(:,pp));
    end;
    fprintf(fid,']'';\n\n');
end;

if exist('imstrnum','var'),
    string_save = [string_save ' imbase imformat imstrnum'];
else
    fprintf(1,'\nNo image found in the calibration results!\n');
end;

if exist('x_cell','var'),
    if ~exist('err_std','var'),
        fprintf(1,'\nThe pixel error not found, you need to calibrate first (press "Calibration button")!\n');
        fclose(fid);
        return;
    end;

    if ~exist('map','var'),
        map = gray(256);
    end;

    if ~exist('est_fc_mat','var');
        est_fc_mat = ones(2,n_cam);    % Set to zero if you do not want to estimate the focal length
    end;

    if ~exist('center_optim_vec','var'),
        center_optim_vec = ones(1,n_cam);
    end;

    if ~exist('est_aspect_ratio_vec','var'),
        est_aspect_ratio_vec = ones(1,n_cam);
    end;

    if ~exist('est_alpha_vec','var'),
        est_alpha_vec = zeros(1,n_cam);
    end;

    if ~exist('est_dist_mat','var'),
        est_dist_mat = zeros(5,n_cam);
    else
        [m,n] = size(est_dist_mat);
        if n~=n_cam,
            fprintf(1,'\nWARNING: Estimation indicator matrix of distortion do not match number of cameras!\n');
            est_dist_mat = zeros(5,n_cam);
        elseif m<5,
            est_dist_mat = [est_dist_mat;zeros(5-m,n_cam)];
        elseif m>5,
            est_dist_mat = est_dist_mat(1:5,:);
        end;
    end;

    if ~exist('check_cond','var'),
        check_cond = 1;
    end;

    if ~exist('MaxIter','var'),
        MaxIter = 30;
    end;

    fprintf(1,['\nSaving extrinsic calibration results under ' save_name '.mat\n']);
    string_save = [string_save ' n_ima map ind_active active_images err_cam err_std dX dY' ...
                               ' active_imgviews wintx winty win_size est_fc_mat center_optim_vec' ...
                               ' est_alpha_vec est_dist_mat est_aspect_ratio_vec check_cond MaxIter' ...
                               ' X_cell x_cell y_cell ex_cell H_cell Tw_mat_error Tw_mat n_sq_mat' ...
                               ' dXY_mat dX_default dY_default'];

    if rotflag,       % quaternion
        string_save = [string_save ' Qw_mat Qw_mat_error'];
    else              % axis angle
        string_save = [string_save ' Omw_mat Omw_mat_error'];
    end;

    fprintf(fid,'%%%% -- Various other variables (may be ignored if you do not use the Matlab Calibration Toolbox):\n');
    fprintf(fid,'%%-- Those variables are used to control which intrinsic parameters should be optimized\n\n');
    fprintf(fid,'%% Number of image frames:\n');
    fprintf(fid,'n_ima = %d;\n\n',n_ima);

    fprintf(fid,'%% Estimation indicator matrix of the focal length:\n');
    fprintf(fid,'est_fc_mat = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%d, %d;\n', est_fc_mat(:,pp));
    end;
    fprintf(fid,']'';\n\n');

    fprintf(fid,'%% Estimation indicator vector of the aspect ratio fc(2)/fc(1):\n');
    fprintf(fid,'est_aspect_ratio_vec = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%d;\n', est_aspect_ratio_vec(pp));
    end;
    fprintf(fid,']'';\n\n');

    fprintf(fid,'%% Estimation indicator vector of the principal point:\n');
    fprintf(fid,'center_optim_vec = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%d;\n', center_optim_vec(pp));
    end;
    fprintf(fid,']'';\n\n');

    fprintf(fid,'%% Estimation indicator vector of the skew coefficient:\n');
    fprintf(fid,'est_alpha_vec = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%d;\n', est_alpha_vec(pp));
    end;
    fprintf(fid,']'';\n\n');

    fprintf(fid,'%% Estimation indicator matrix of the distortion coefficients:\n');
    fprintf(fid,'est_dist_mat = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%d, %d, %d, %d, %d;\n',est_dist_mat(:,pp));
    end;
    fprintf(fid,']'';\n\n\n');

    if n_cam>1,
        if ~exist('refine_multicam','var'),
            refine_multicam = 1;
        end;
        string_save = [string_save ' refine_multicam'];

        if refine_multicam,
            fprintf(1,['\nWith ' num2str(n_cam) ' cameras available, Calibration results are refined!\n']);
            fprintf(fid,'%%%% -- Refined extrinsic parameters of camera 1 wrt world reference frame:\n');

            if rotflag,
                fprintf(fid,'%%-- Quaternion rotation (Qcw) and translation (Tcw) vectors and their uncertainties\n\n');
                string_save = [string_save ' extrinsic_deviation Qcw Tcw Qcw_error Tcw_error'];
                fprintf(fid,'Qcw = [\n');
                for kk = 1:n_ima,
                    fprintf(fid,'%5.15f, %5.15f, %5.15f, %5.15f;\n', Qcw(:,kk));
                end;
                fprintf(fid,']'';\n\n');

                fprintf(fid,'Qcw_error = [\n');
                for kk = 1:n_ima,
                    fprintf(fid,'%5.15f, %5.15f, %5.15f, %5.15f;\n', Qcw_error(:,kk));
                end;
                fprintf(fid,']'';\n\n');
            else
                fprintf(fid,'%%-- Axis angle rotation (Omcw) and translation (Tcw) vectors and their uncertainties\n\n');
                string_save = [string_save ' extrinsic_deviation Omcw Tcw Omcw_error Tcw_error'];
                fprintf(fid,'Omcw = [\n');
                for kk = 1:n_ima,
                    fprintf(fid,'%5.15f, %5.15f, %5.15f;\n', Omcw(:,kk));
                end;
                fprintf(fid,']'';\n\n');

                fprintf(fid,'Omcw_error = [\n');
                for kk = 1:n_ima,
                    fprintf(fid,'%5.15f, %5.15f, %5.15f;\n', Omcw_error(:,kk));
                end;
                fprintf(fid,']'';\n\n');
            end;

            fprintf(fid,'Tcw = [\n');
            for kk = 1:n_ima,
                fprintf(fid,'%5.15f, %5.15f, %5.15f;\n', Tcw(:,kk));
            end;
            fprintf(fid,']'';\n\n');

            fprintf(fid,'Tcw_error = [\n');
            for kk = 1:n_ima,
                fprintf(fid,'%5.15f, %5.15f, %5.15f;\n', Tcw_error(:,kk));
            end;
            fprintf(fid,']'';\n\n\n');
        end;
    end;

    fprintf(fid,'%%%% -- Extrinsic parameters of each camera wrt world reference frame:\n');
    if ~exist('n_view','var'),
        n_view = n_ima * n_cam;
    end;

    if rotflag,
        fprintf(fid,'%%-- Quaternion rotation (Qw_mat) and translation (Tw_mat) vectors and their uncertainties\n\n');
        fprintf(fid,'Qw_mat = [\n');
        for kk = 1:n_view,
            fprintf(fid,'%5.15f, %5.15f, %5.15f, %5.15f;\n', Qw_mat(:,kk));
        end;
        fprintf(fid,']'';\n\n');

        fprintf(fid,'Qw_mat_error = [\n');
        for kk = 1:n_view,
            fprintf(fid,'%5.15f, %5.15f, %5.15f, %5.15f;\n', Qw_mat_error(:,kk));
        end;
        fprintf(fid,']'';\n\n');
    else
        fprintf(fid,'%%-- Axis angle rotation (Omw_mat) and translation (Tw_mat) vectors and their uncertainties\n\n');
        fprintf(fid,'Omw_mat = [\n');
        for kk = 1:n_view,
            fprintf(fid,'%5.15f, %5.15f, %5.15f;\n', Omw_mat(:,kk));
        end;
        fprintf(fid,']'';\n\n');

        fprintf(fid,'Omw_mat_error = [\n');
        for kk = 1:n_view,
            fprintf(fid,'%5.15f, %5.15f, %5.15f;\n', Omw_mat_error(:,kk));
        end;
        fprintf(fid,']'';\n\n');
    end;

    fprintf(fid,'Tw_mat = [\n');
    for kk = 1:n_view,
        fprintf(fid,'%5.15f, %5.15f, %5.15f;\n', Tw_mat(:,kk));
    end;
    fprintf(fid,']'';\n\n');

    fprintf(fid,'Tw_mat_error = [\n');
    for kk = 1:n_view,
        fprintf(fid,'%5.15f, %5.15f, %5.15f;\n', Tw_mat_error(:,kk));
    end;
    fprintf(fid,']'';\n\n');
end;

eval(string_save);
fclose(fid);
fprintf(1,'done...To load later click on Load\n');
