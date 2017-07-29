%function multicams_gui,
cell_list = cell(4,4);
if ~exist('rotflag','var'),
    rotflag = 0;
end;

%-------- Begin editable region -------------%
fig_number = 1;
title_figure = 'Multicameras Calibration Toolbox --- Different Intrinsic';

cell_list{1,1} = {'Read Images','read_images_multicam;'};
cell_list{1,2} = {'Extract Grid Corners','click_track_multicam;'};
cell_list{1,3} = {'Recheck Corners','recheck_corners_multicam;'};
if rotflag,
    cell_list{1,4} = {'Calibration','go_calib_optim_multicam;'};       % quaternion rotation (need correction)
else
    cell_list{1,4} = {'Calibration','go_calib_optim_multicam2;'};     % axis angle rotation
end;
cell_list{2,1} = {'Convert Format','convert_bouguet_multicam;'};
cell_list{2,2} = {'Analyse Error','analyse_error_multicam;'};
cell_list{2,3} = {'Reproject on Images','reproject_calib_multicam;'};
cell_list{2,4} = {'Recompute Corners','recompute_corners_multicam;'};
cell_list{3,1} = {'Show Intrinsic','show_intrinsic_multicam;'};
cell_list{3,2} = {'Show Extrinsic','plot_camera_multicam;'};
cell_list{3,3} = {'Visualize Distortion','visualize_distortion_multicam;'};
cell_list{3,4} = {'Undistort Images','undistort_image_multicam;'};
cell_list{4,1} = {'Save','saving_calib_multicam;'};
cell_list{4,2} = {'Load','loading_calib_multicam;'};
cell_list{4,3} = {'Export Data','export_calib_multicam;'};
cell_list{4,4} ={'Exit',['disp(''Bye. To run again, type multiview_gui.''); close(' num2str(fig_number) ');']};

show_window(cell_list,fig_number,title_figure,135,18,0,'clean',12);
