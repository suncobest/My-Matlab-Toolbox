%function multiview_gui(mode)

% multiview_gui(mode)
%
% Runs the extended Camera Calibration Toolbox for multiple views.
%
% TWO MODES OF THE MULTIVIEW CALIBRATION TOOLBOX:
% Two versions of Multiview application: 1. Multiple cameras; 2. One amera with split views.
clc;
cell_list = {4,3};
fig_number = 1;

title_figure = 'Visual Hull Construction Toolbox:';
cell_list{1,1} = {'Load cameras','load_cameras;'};
cell_list{1,2} = {'Foreground extraction','foreground_extraction;'};
cell_list{1,3} = {'3D reconstruction','visualhull_carving;'};
cell_list{2,1} = {'Show structure','show_structure;'};
cell_list{2,2} = {'Segment structure','segment_structure;'};
cell_list{2,3} = {'Refine segmentation','refine_segment;'};
cell_list{3,1} = {'Kinematics extraction','kinematics_extraction;'};
cell_list{3,2} = {'Polygon animation','polygon_animation;'};
cell_list{3,3} = {'Check reprojection','check_reprojection;'};

% set global ColorOrder
palette = [0, 0, 1; 0, 0.5, 0; 1, 0, 0; 0, 0.75, 0.75; 0.75, 0, 0.75; 0.75, 0.75, 0; 0.25, 0.25, 0.25];
set(0,'DefaultAxesColorOrder',palette);

show_window(cell_list,fig_number,title_figure,150,18,0,'clean',12);