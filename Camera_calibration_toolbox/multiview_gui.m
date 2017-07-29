%function multiview_gui(mode)

% multiview_gui(mode)
%
% Runs the extended Camera Calibration Toolbox for multiple views.
%
% TWO MODES OF THE MULTIVIEW CALIBRATION TOOLBOX:
% Two versions of Multiview application: 1. Multiple cameras; 2. Camera with mirros.

clear;
cell_list = {};
fig_number = 1;

% fprintf(1,'Which kind of rotation vector do you want to use?\n');
% rotflag = input('Rotation vector: ([]=Axis angle, other=Quaternion) ','s');
% rotflag = ~isempty(rotflag);

title_figure = 'Multiview Calibration Toolbox --- Select Mode:';

cell_list{1,1} = {'Multiple Cameras Configuration','multicams_gui;'};
cell_list{2,1} = {'One Camera with Split Views','multimirrors_gui;'};
cell_list{3,1} = {'Exit',['disp(''Bye. To run again, type multiview_gui.''); close(' num2str(fig_number) ');']};


show_window(cell_list,fig_number,title_figure,250,20,0,'clean',12);