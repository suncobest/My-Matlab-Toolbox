%function Plane_calib

cell_list = {};

%-------- Begin editable region -------------%

fig_number = 1;
title_figure = 'Calibrate Plane ';


cell_list{1,1} = {'Read Image','read_image;'};
cell_list{2,1} = {'Calibrate With Triangle','calibrate_plane_tri;'};
cell_list{3,1} = {'Calibrate With Rectangle','calibrate_plane_rect;'};
cell_list{4,1} = {'Save Calibration Result','save_calibration_result;'};
cell_list{5,1} = {'Load Calibration Result','load_calibration_result;'};
cell_list{6,1} = {'Measure Length','measure_length;'};
cell_list{7,1} = {'Measure Angle','measure_angle;'};


show_window(cell_list,fig_number,title_figure,190,20,0,'clean',12);


%-------- End editable region -------------%


