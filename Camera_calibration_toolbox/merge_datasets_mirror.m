% call_flag is set in case the script called by loading_calib_mirror
if ~exist('call_flag','var'),
    call_flag = 0;
end;
if ~call_flag,
    fprintf(1,'Script that merges multiple "Multimirror_Calib_Results*.mat" datasets of same cameras into a single dataset.\n');
    dir('Multimirror_Calib_Results*.mat');
    datasets = dir('Multimirror_Calib_Results*.mat');
    Ndata = length(datasets);  % number of Multimirror_Calib_Results*.mat
    if Ndata==0,
        fprintf(1,'Multimirror_Calib_Results*.mat not found in current directory!\n');
        return;
    elseif Ndata==1,
        fprintf(1,'Only one Multimirror_Calib_Results*.mat file found in current directory! Load file ...\n');
        load(datasets(1).name);
        return;
    end;
end;

fprintf(1,'Make sure all "Multimirror_Calib_Results*.mat" datasets under current path from same camera system!\n');
% load Multimirror_Calib_Results
% common variables
dataii = load(datasets(1).name);
nx = dataii.nx;
ny = dataii.ny;
n_cam = dataii.n_cam;
hand_list = dataii.hand_list;

if isfield(dataii,'fc'),
    fc = dataii.fc;
    cc = dataii.cc;
    kc = dataii.kc;
    alpha_c = dataii.alpha_c;
end;

fprintf(1,'\nThe estimation indicators of merging data will be same as %s ...\n',datasets(1).name);
if isfield(dataii,'est_fc'),
    est_fc = dataii.est_fc;
end;

if isfield(dataii,'center_optim'),
    center_optim = dataii.center_optim;
end;

if isfield(dataii,'est_dist'),
    est_dist = dataii.est_dist;
end;

if isfield(dataii,'est_alpha'),
    est_alpha = dataii.est_alpha;
end;

if isfield(dataii,'est_aspect_ratio'),
    est_aspect_ratio = dataii.est_aspect_ratio;
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
if isfield(dataii,'strnum_cell'),
    calib_name = dataii.calib_name;
    format_image = dataii.format_image;
    flag = 1;
end;

% check data
n_ima = dataii.n_ima;
for ii = 2:Ndata,
    dataii = load(datasets(ii).name);
    if flag,
        if ~strcmp(calib_name, dataii.calib_name),
            fprintf(1,'\nWARMING: Image basename in data #%d do not match with data #1!\n',ii);
            flag = 0;
        end;
        if ~strcmp(format_image, dataii.format_image),
            fprintf(1,'\nWARMING: Image suffix in data #%d do not match with data #1!\n',ii);
            flag = 0;
        end;
    end;
    if nx ~= dataii.nx || ny ~= dataii.ny,
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
    strnum_cell = cell(1,n_ima);
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
        strnum_cell(indii+1 : indii+nimii) = dataii.strnum_cell;
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
fprintf('Multi-mirror datasets are now merged. You are now ready to run calibration. \n');

% datasets = rmfield(datasets,'date');
% datasets = rmfield(datasets,'bytes');
% datasets = rmfield(datasets,'isdir');
% datasets = rmfield(datasets,'datenum');
