% call_flag is set in case the script called by loading_calib_multicam
if ~exist('call_flag','var'),
    call_flag = 0;
end;
if ~call_flag,
    fprintf(1,'Script that merges multiple "Multicam_Calib_Results*.mat" datasets of same cameras into a single dataset.\n');
    dir('Multicam_Calib_Results*.mat');
    datasets = dir('Multicam_Calib_Results*.mat');
    Ndata = length(datasets);  % number of Multicam_Calib_Results*.mat
    if Ndata==0,
        fprintf(1,'Multicam_Calib_Results*.mat not found in current directory!\n');
        return;
    elseif Ndata==1,
        fprintf(1,'Only one Multicam_Calib_Results*.mat file found in current directory! Load file ...\n');
        load(datasets(1).name);
        return;
    end;
end;

fprintf(1,'Make sure all "Multicam_Calib_Results*.mat" datasets under current path from same camera system!\n');
% load Multicam_Calib_Results
% common variables
dataii = load(datasets(1).name);
imsize = dataii.imsize;
n_cam = dataii.n_cam;
hand_list = dataii.hand_list;

if isfield(dataii,'fc_mat'),
    fc_mat = dataii.fc_mat;
    cc_mat = dataii.cc_mat;
    kc_mat = dataii.kc_mat;
    alpha_vec = dataii.alpha_vec;
end;

fprintf(1,'\nThe estimation indicators of merging data will be same as %s ...\n',datasets(1).name);
if isfield(dataii,'est_fc_mat'),
    est_fc_mat = dataii.est_fc_mat;
end;

if isfield(dataii,'center_optim_vec'),
    center_optim_vec = dataii.center_optim_vec;
end;

if isfield(dataii,'est_dist_mat'),
    est_dist_mat = dataii.est_dist_mat;
end;

if isfield(dataii,'est_alpha_vec'),
    est_alpha_vec = dataii.est_alpha_vec;
end;

if isfield(dataii,'est_aspect_ratio_vec'),
    est_aspect_ratio_vec = dataii.est_aspect_ratio_vec;
end;

if isfield(dataii,'check_cond'),
    check_cond = dataii.check_cond;
end;

if isfield(dataii,'MaxIter'),
    MaxIter = dataii.MaxIter;
end;

if isfield(dataii,'map'),
    map = dataii.map;
end;

if isfield(dataii,'dX_default'),
    dX = dataii.dX;
    dY = dataii.dY;
    dX_default = dataii.dX_default;
    dY_default = dataii.dY_default;
end;

if isfield(dataii,'win_size'),
    wintx = dataii.wintx;
    winty = dataii.winty;
    win_size = dataii.win_size;
end;

flag = 0;
if isfield(dataii,'imstrnum'),
    imbase = dataii.imbase;
    imformat = dataii.imformat; 
    flag = 1;
end;

% check data
n_ima = dataii.n_ima;
for ii = 2:Ndata,
    dataii = load(datasets(ii).name);
    if flag,
        if ~isequal(imbase, dataii.imbase),
            fprintf(1,'\nWARMING: Image basename in data #%d do not match with data #1!\n',ii);
            flag = 0;
        end;
        if ~isequal(imformat, dataii.imformat),
            fprintf(1,'\nWARMING: Image suffix in data #%d do not match with data #1!\n',ii);
            flag = 0;
        end;
    end;
    if any(imsize(:) ~= dataii.imsize(:)),
        fprintf(1,'\nERROR: Picture sizes in your data do not match with each other!\n');
        return;
    end;
    if n_cam ~= dataii.n_cam,
        fprintf(1,'\nERROR: Number of cameras in your data do not match with each other!\n');
        return;
    end;
    if any(hand_list ~= dataii.hand_list),
        fprintf(1,'\nERROR: Handedness of your cameras do not match with each other!\n');
        return;
    end;
    n_ima = n_ima + dataii.n_ima;
end;

if flag,
    imstrnum = cell(1,n_cam);
    for pp = 1:n_cam,
        imstrnum{pp} = cell(1,n_ima);
    end;
end;
% variables to update
n_view = n_ima * n_cam;
active_imgviews = zeros(n_cam,n_ima);
X_cell =  cell(1,n_view);
x_cell =  X_cell;
dXY_mat = NaN(2,n_view);
n_sq_mat =  dXY_mat;
indii = 0;
for ii = 1:Ndata,
    dataii = load(datasets(ii).name);
    nimii = dataii.n_ima;
    if flag,
        for pp=1:n_cam,
            imstrnum{pp}(indii+1 : indii+nimii) = dataii.imstrnum{pp};
        end;
    end;
    active_imgviews(:,indii+1 : indii+nimii) = dataii.active_imgviews;
    kth = indii*n_cam+1 : (indii+nimii)*n_cam;
    X_cell(kth) =  dataii.X_cell;
    x_cell(kth) = dataii.x_cell;
    dXY_mat(:,kth) = dataii.dXY_mat;
    n_sq_mat(:,kth) = dataii.n_sq_mat;
    indii = indii+nimii;
end;
active_images = any(active_imgviews,1);
ind_active = find(active_images);

clear dataii datasets
fprintf('Multi-camera datasets are now merged. You are now ready to run calibration. \n');