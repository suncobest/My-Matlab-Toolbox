if ~exist('n_ima','var')||~exist('fc_mat','var')||~exist('cc_mat','var'),
    fprintf(1,'No calibration data available.\n');
    return;
end;

if ~exist('n_cam','var'),
    n_cam = 1;
end;

if ~exist('alpha_vec','var'),
    alpha_vec = zeros(1,n_cam);
end;

if ~exist('kc_mat','var'),
    kc_mat = zeros(5,n_cam);
end;

if ~exist('fc_mat_error','var'),
    for pp = 1:n_cam,
        fc = fc_mat(:,pp);
        cc = cc_mat(:,pp);
        kc = kc_mat(:,pp);
        alpha_c = alpha_vec(pp);
        fprintf(1,'\n\nCalibration results of camera %d:\n\n',pp);
        fprintf(1,'Focal Length:          fc = [%3.5f, %3.5f]\n',fc);
        fprintf(1,'Principal point:       cc = [%3.5f, %3.5f]\n',cc);
        fprintf(1,'Skew:             alpha_c =  %3.5f  => Skew angle = %3.5f degrees\n',alpha_c,90-atan(alpha_c)*180/pi);
        fprintf(1,'Distortion:              kc = [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f]\n',kc);
    end;
else
    for pp = 1:n_cam,
        fc = fc_mat(:,pp);
        cc = cc_mat(:,pp);
        kc = kc_mat(:,pp);
        alpha_c = alpha_vec(pp);
        fc_error = fc_mat_error(:,pp);
        cc_error = cc_mat_error(:,pp);
        kc_error = kc_mat_error(:,pp);
        alpha_c_error = alpha_vec_error(pp);
        fprintf(1,'\n\nCalibration results of camera %d (with uncertainties):\n\n',pp);
        fprintf(1,'Focal Length:      fc = [%3.5f, %3.5f] ? [%3.5f, %3.5f]\n',[fc;fc_error]);
        fprintf(1,'Principal point:   cc = [%3.5f, %3.5f] ? [%3.5f, %3.5f]\n',[cc;cc_error]);
        fprintf(1,'Skew:         alpha_c = %3.5f ? %3.5f => Skew angle = %3.5f ? %3.5f degrees\n', ...
            [alpha_c; alpha_c_error; 90-atan(alpha_c)*180/pi; atan(alpha_c_error)*180/pi]);
        fprintf(1,'Distortion:          kc = [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f] ? [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f]\n',[kc;kc_error]);
        if n_ima ~= 0,
            fprintf(1,'Pixel error:        err = [%3.5f, %3.5f]\n\n',err_cam(:,pp));
        end;
        fprintf(1,['\nNote: The numerical errors are approximately three times ' ...
            'the standard deviations (for reference).\n\n']);
    end;
end;
fprintf(1,'Done!\n\n');
