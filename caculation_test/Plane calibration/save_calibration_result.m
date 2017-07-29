if ~exist('mag','var')
   fprintf(1,'No calibration data available.\n');
   return;
end;

save calibrated_plane.mat I wintx winty focus pix mag dir_plane theta so si
if exist('Homo','var')
    save calibrated_plane.mat Homo
end

if exist('oxy','var')
    save calibrated_plane.mat oxy
end

disp('calibration data saved!')