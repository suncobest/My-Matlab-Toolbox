if ~exist('n_ima','var')||~exist('fc','var')||~exist('cc','var'),
    fprintf(1,'No calibration data available.\n');
    return;
end;

if ~exist('alpha_c','var'),
    alpha_c = 0;
end;

if ~exist('kc','var'),
    kc = zeros(5,1);
elseif length(kc) == 4;
    kc = [kc;0];
end;

if ~exist('fc_error','var'),
    fprintf(1,'\n\nCalibration results:\n\n');
    fprintf(1,'Focal Length:          fc = [%3.5f, %3.5f]\n',fc);
    fprintf(1,'Principal point:       cc = [%3.5f, %3.5f]\n',cc);
    fprintf(1,'Skew:            alpha_c =  %3.5f  => Skew angle = %3.5f degrees\n',alpha_c,90-atan(alpha_c)*180/pi);
    fprintf(1,'Distortion:             kc = [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f]\n',kc);
else
    fprintf(1,'\n\nCalibration results (with uncertainties):\n\n');
    fprintf(1,'Focal Length:       fc = [%3.5f, %3.5f] ? [%3.5f, %3.5f]\n',[fc;fc_error]);
    fprintf(1,'Principal point:     cc = [%3.5f, %3.5f] ? [%3.5f, %3.5f]\n',[cc;cc_error]);
    fprintf(1,'Skew:         alpha_c = %3.5f ? %3.5f => Skew angle = %3.5f ? %3.5f degrees\n', ...
        [alpha_c; alpha_c_error; 90-atan(alpha_c)*180/pi; atan(alpha_c_error)*180/pi]);
    fprintf(1,'Distortion:           kc = [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f] ? [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f]\n',[kc;kc_error]);
    if n_ima ~= 0,
        fprintf(1,'Pixel error of every camera view:\n');
        for pp = 1:n_cam,
            fprintf(1,'   view %d:         err = [%3.5f, %3.5f]\n',pp,err_cam(:,pp));
        end;
    end;
    fprintf(1,['\nNote: The numerical errors are approximately three times ' ...
        'the standard deviations (for reference).\n\n']);
end;
fprintf(1,'Done!\n\n');
