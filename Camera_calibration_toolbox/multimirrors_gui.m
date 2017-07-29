%function multimirrors_gui,
cell_list = cell(4,4);
if ~exist('rotflag','var'),
    rotflag = 0;
end;

%-------- Begin editable region -------------%
fig_number = 1;
title_figure = 'Multimirrors Calibration Toolbox --- Identical Intrinsic';

cell_list{1,1} = {'Read Images','read_images_mirror;'};
cell_list{1,2} = {'Extract Grid Corners','click_track_mirror;'};
cell_list{1,3} = {'Recheck Corners','recheck_corners_mirror;'};
if rotflag,
    cell_list{1,4} = {'Calibration','go_calib_optim_mirror;'};       % quaternion rotation (need correction)
else
    cell_list{1,4} = {'Calibration','go_calib_optim_mirror2;'};     % axis angle rotation
end;
cell_list{2,1} = {'Convert Format','convert_bouguet_mirror;'};
cell_list{2,2} = {'Analyse Error','analyse_error_mirror;'};
cell_list{2,3} = {'Reproject on Images','reproject_calib_mirror;'};
cell_list{2,4} = {'Recompute Corners','recompute_corners_mirror;'};
cell_list{3,1} = {'Show Intrinsic','show_intrinsic_mirror;'};
cell_list{3,2} = {'Show Extrinsic','plot_camera_mirror;'};
cell_list{3,3} = {'Visualize Distortion','visualize_distortion_mirror;'};
cell_list{3,4} = {'Undistort Images','undistort_image_mirror;'};
cell_list{4,1} = {'Save','saving_calib_mirror;'};
cell_list{4,2} = {'Load','loading_calib_mirror;'};
cell_list{4,3} = {'Export Data','export_calib_mirror;'};
cell_list{4,4} ={'Exit',['disp(''Bye. To run again, type multiview_gui.''); close(' num2str(fig_number) ');']};

show_window(cell_list,fig_number,title_figure,135,18,0,'clean',12);
