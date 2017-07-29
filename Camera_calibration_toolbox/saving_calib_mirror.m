if ~exist('cc','var')||~exist('fc','var'),
    fprintf(1,'No calibration data available.\n');
    return;
end;

if ~exist('n_cam','var'),
    n_cam = 1;
end;

if ~exist('kc','var'),
    kc = zeros(5,1);
else
    kc = kc(:);
    n = length(kc);
    if n<5,
        kc = [kc;zeros(5-n,1)];
    elseif n>5,
        kc = kc(1:5);
    end;
end;

if ~exist('alpha_c','var'),
    alpha_c = 0;
end;

if ~exist('hand_list','var'),
    hand_list = ones(1,n_cam);
end;

save_name = 'Multimirror_Calib_Results';

if exist([save_name '.mat'],'file')==2,
    disp('WARNING: File ''Multimirror_Calib_Results.mat'' already exists!');
    if exist('copyfile','builtin')==5,
        cont = -1;
        flag = 1;
        while flag,
            cont = cont + 1;
            save_name = ['Old_Multimirror_Calib_Results' num2str(cont)];
            flag = (exist([save_name '.mat'],'file')==2);
        end;
        copyfile('Multimirror_Calib_Results.mat',[save_name '.mat']);
        disp(['Copying the current ''Multimirror_Calib_Results.mat'' file to ''' save_name '.mat''']);
        if exist('Multimirror_Calib_Results.m','file')==2,
            copyfile('Multimirror_Calib_Results.m',[save_name '.m']);
            disp(['Copying the current ''Multimirror_Calib_Results.m'' file to ''' save_name '.m''']);
        end;
        cont_save = 1;
    else
        disp('The file ''Multimirror_Calib_Results.mat'' is about to be changed.');
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

save_name = 'Multimirror_Calib_Results';
fprintf(1,['\nSaving intrinsic calibration results under ''' save_name '.mat''\n']);
string_save = ['save ' save_name ' fc cc kc alpha_c n_cam hand_list nx ny'];

fid = fopen([save_name '.m'],'wt');
fprintf(1,'Generating the matlab script file ''%s.m'' containing the intrinsic and extrinsic parameters...\n',save_name);
fprintf(fid,'%%%% Intrinsic and Extrinsic Camera Parameters:\n');
fprintf(fid,'%%\n');
fprintf(fid,'%% This script file can be directly excecuted under Matlab to recover the camera intrinsic and extrinsic parameters.\n');
fprintf(fid,'%% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.\n');
fprintf(fid,'%%            All those complementary variables are saved in the complete matlab data file Multimirror_Calib_Results.mat.\n');
fprintf(fid,'%% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/doc/\n\n\n');

fprintf(fid,'%%-- Focal length:\n');
fprintf(fid,'fc = [ %5.15f ; %5.15f ];\n\n',fc);

fprintf(fid,'%%-- Principal point:\n');
fprintf(fid,'cc = [ %5.15f ; %5.15f ];\n\n',cc);

fprintf(fid,'%%-- Skew coefficient:\n');
fprintf(fid,'alpha_c = %5.15f;\n\n',alpha_c);

fprintf(fid,'%%-- Distortion coefficients:\n');
fprintf(fid,'kc = [ %5.15f ; %5.15f ; %5.15f ; %5.15f ; %5.15f ];\n\n',kc);

fprintf(fid,'%%-- Number of views in every image:\n');
fprintf(fid,'n_cam = %d;\n\n',n_cam);

fprintf(fid,'%%-- Axis handness of all camera views:\n');
fprintf(fid,'hand_list = [\n');
for pp = 1:n_cam,
    fprintf(fid,'%d;\n', hand_list(pp));
end;
fprintf(fid,']'';\n\n');

fprintf(fid,'%%-- Image size:\n');
fprintf(fid,'nx = %d;\n',nx);
fprintf(fid,'ny = %d;\n\n',ny);

if exist('fc_error','var'),
    fprintf(1,['\nSaving intrinsic error under ' save_name '.mat\n']);
    string_save = [string_save ' fc_error cc_error alpha_c_error kc_error'];

    fprintf(fid,'%%-- Focal length uncertainty:\n');
    fprintf(fid,'fc_error = [ %5.15f ; %5.15f ];\n\n',fc_error);

    fprintf(fid,'%%-- Principal point uncertainty:\n');
    fprintf(fid,'cc_error = [ %5.15f ; %5.15f ];\n\n',cc_error);

    fprintf(fid,'%%-- Skew coefficient uncertainty:\n');
    fprintf(fid,'alpha_c_error = %5.15f;\n\n',alpha_c_error);

    fprintf(fid,'%%-- Distortion coefficients uncertainty:\n');
    fprintf(fid,'kc_error = [ %5.15f ; %5.15f ; %5.15f ; %5.15f ; %5.15f ];\n\n\n',kc_error);
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
    fprintf(1,'\nConfiguration of multiple views detected!\n');

    fprintf(fid,'%%-- Axis handness of all camera views wrt view 1:\n');
    fprintf(fid,'handcc = [\n');
    for pp = 1:n_cam,
        fprintf(fid,'%d;\n', handcc(pp));
    end;
    fprintf(fid,']'';\n\n');

    fprintf(fid,'%%-- Refined extrinsic parameters of all camera views wrt view 1:\n');
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

if exist('strnum_cell','var'),
    string_save = [string_save ' calib_name format_image strnum_cell'];
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

    if ~exist('est_fc','var');
        est_fc = [1;1]; % Set to zero if you do not want to estimate the focal length
    end;

    if ~exist('center_optim','var'),
        center_optim = 1;
    end;

    if ~exist('est_aspect_ratio','var'),
        est_aspect_ratio = 1;
    end;

    if ~exist('est_alpha','var'),
        est_alpha = 0;
    end;

    if ~exist('est_dist','var'),
        est_dist = zeros(5,1);
    else
        n = length(est_dist);
        if n<5,
            kc = [est_dist;zeros(5-n,1)];
        elseif n>5,
            est_dist = est_dist(1:5);
        end;
    end;

    if ~exist('check_cond','var'),
        check_cond = 1;
    end;

    if ~exist('MaxIter','var'),
        MaxIter = 30;
    end;

    fprintf(1,['\nSaving extrinsic calibration results under ' save_name '.mat\n']);
    string_save = [string_save ' n_ima map ind_active active_images err_cam err_std' ...
                               ' active_imgviews wintx winty win_size est_fc center_optim' ...
                               ' est_alpha est_dist est_aspect_ratio check_cond MaxIter' ...
                               ' X_cell x_cell y_cell ex_cell H_cell Tw_mat_error Tw_mat' ...
                               ' n_sq_mat dXY_mat dX_default dY_default dX dY'];

    if rotflag,     % quaternion
        string_save = [string_save ' Qw_mat Qw_mat_error'];
    else             % axis angle
        string_save = [string_save ' Omw_mat Omw_mat_error'];
    end;

    fprintf(fid,'%%%% -- Various other variables (may be ignored if you do not use the Matlab Calibration Toolbox):\n');
    fprintf(fid,'%%-- Those variables are used to control which intrinsic parameters should be optimized\n\n');
    fprintf(fid,'%% Number of image frames:\n');
    fprintf(fid,'n_ima = %d;\n',n_ima);
    fprintf(fid,'%% Estimation indicator of the two focal variables:\n');
    fprintf(fid,'est_fc = [ %d ; %d ];\n',est_fc);
    fprintf(fid,'%% Estimation indicator of the aspect ratio fc(2)/fc(1):\n');
    fprintf(fid,'est_aspect_ratio = %d;\n',est_aspect_ratio);
    fprintf(fid,'%% Estimation indicator of the principal point:\n');
    fprintf(fid,'center_optim = %d;\n',center_optim);
    fprintf(fid,'%% Estimation indicator of the skew coefficient:\n');
    fprintf(fid,'est_alpha = %d;\n',est_alpha);
    fprintf(fid,'%% Estimation indicator of the distortion coefficients:\n');
    fprintf(fid,'est_dist = [ %d; %d ; %d ; %d ; %d ];\n\n\n',est_dist);

    if n_cam>1,
        if ~exist('refine_multicam','var'),
            refine_multicam = 1;
        end;
        string_save = [string_save ' refine_multicam'];

        if refine_multicam,
            fprintf(1,['\nWith ' num2str(n_cam) ' camera views available, Calibration results are refined!\n']);
            fprintf(fid,'%%%% -- Refined extrinsic parameters of view 1 wrt world reference frame:\n');

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

    fprintf(fid,'%%%% -- Extrinsic parameters of each view wrt world reference frame:\n');
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
