if exist('calibrated_plane.mat','file')
    load('calibrated_plane.mat')
    disp('Calibration result (calibrated_plane.mat) loaded!')
else
    error('No calibration result (calibrated_plane.mat) can be found!')
end