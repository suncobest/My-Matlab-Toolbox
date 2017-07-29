% find the codedir of the toolbox
codeDir = mfilename('fullpath');
flag = find(codeDir=='\' | codeDir=='/',1,'last');
if ~isempty(flag),
    codeDir = codeDir(1 : flag(end)-1);
end;
if exist([codeDir '/load_imgdir.m'],'file'),
    load_imgdir;
    fprintf(1,'\nModify file ''load_imgdir.m'' in code directory if you want to redirect image path!\n');
else
    fprintf(1,'\nNo file ''load_imgdir.m'' found in the code directory!\nPlease compute visual hull first!\n');
    return;
end;

save_name = [imgdir '/kinematic_animation.mat'];
if exist(save_name,'file')==2,
    fprintf(1,'\nLoading kinematics parameters file ''kinematic_animation.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\nERROR: ''kinematic_animation.mat'' not found!\n');
    return;
end;

save_name = [imgdir '/visualhull_environment.mat'];
if exist(save_name, 'file')==2,
    fprintf(1,'\nLoading environment variable file ''visualhull_environment.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\nERROR: ''visualhull_environment.mat'' not found!\n');
    return;
end;

% initial 2*2 matrix for perspective image corner (for surface function)
xc = [0 1;0 1];
yc = [0 0;1 1];
zc = [1 1;1 1];

% view axis for 3D space
az = 50;
el = 30;
palette = lines(7);        % color of faces
palette2 = palette;
% palette2 = char(ones(7,1)*'none');      % color of wire frame

w1 = 0.4; % linewidth of axis in 2D projection
w2 = 0.1; % linewidth of polygon edge in 2D projection
w3 = 1.5; % linewidth of axis in 3D view
al = 0;     % face alpha in 2D prjection (0~1)

if ~exist('map','var'),
    map = gray(256);
end;

nlen_vec = input('Length factor of principle axis vector to show: ([]=1.5 times semi_axes) ');
if isempty(nlen_vec),
    nlen_vec = 1.5;
end;

a = min(member_center(:,1,ind_active),[],3)-maxbound/2;
b = max(member_center(:,1,ind_active),[],3)+maxbound/2;
center = (a+b)/2;     % center of 3D space
ax = max(b-a)/2;
temp = max(ax, img_position);
img_position = temp-maxbound/10;
axes_range = [center(1)-temp,center(1)+temp,center(3)-temp,center(3)+temp,-center(2)-temp,-center(2)+temp];
Xim = zeros(n_cam*3,4);
xim = zeros(n_cam*2,4);
for pp = 1:n_cam,
    om = Omcc(:,pp);
    T = Tcc(:,pp);
    hand = handcc(pp);
    if ~calib_mode,    % multiple cameras
        fc = fc_mat(:,pp);
        cc = cc_mat(:,pp);
        alpha_c = alpha_vec(pp);
    end;
    % direction of camera frame axes in world frame: Rkk
    % camera center in world frame: -Rkk*T
    Rkk = rodrigues(-om);       % transpose
    if hand~=1,
        Rkk(3,:) = -Rkk(3,:);
    end;
    % positon of pespective image center
    vect = center+Rkk*T;
    XX = center+img_position*vect/norm(vect);
    % position of pespective image corners (o, y, xy, x)
    XX = XX(:,ones(1,4))+(Rkk(:,1)*[-1,-1,1,1]+Rkk(:,2)*[-1,1,1,-1])*ax;
    % project the four corner on image
    xx = project_points_mirror2(XX,om,T,hand,fc,cc,zeros(5,1),alpha_c);
    % transform to matlab 3D axes for surface function (o, y, x, xy)
    Xim((pp-1)*3+1 : pp*3, :) = XX(:, [1,2,4,3]);
    xim((pp-1)*2+1 : pp*2, :) = xx;
end;

flag = input('Reset the axes range of 3D scene or not? ([]=no, other=yes) ','s');
set_range = ~isempty(flag);

flag = input('Load STL polygon or ellipsoid to animate insect? ([]=STL model, other=ellipsoid) ','s');
poly_sw = isempty(flag);

% reprojection
if poly_sw,
    func_flag = input('Which function to call: ([]=polygon_projection, other=polygon_center_projection) ','s');
    func_flag = isempty(func_flag);
        % check and load stl model
    load_stl = 1;
    if ~exist('XYZb', 'var') || ~exist('XYZw','var') || isempty(XYZb) || isempty(XYZw),
        load_stl = 0;
    else
        fprintf(1,'\nPolygon variable of body and wing detected!\n');
        flag = input('Do you want to change the model files or not? ([]=no, other=yes) ','s');
        if ~isempty(flag),
            fprintf(1,'\nYou choose to load another STL model ...\n');
            load_stl = 0;
        end;
    end;

    if ~load_stl,
        fprintf(1,'\nPlease make sure the STL models of body and wing are in the directory of ''stlTools''!\n');
        flag = input('Have you prepared STL models yet? ([]=yes, other=no) ','s');
        if ~isempty(flag),
            fprintf(1,'\nSTL models of body and wing are needed afterward...\n');
            return;
        end;
        % find the codedir of the stlTools toolbox
        codeDir = which('stlRead');
        if isempty(codeDir),
            fprintf(1,'\nYou must add the toolbox ''stlTools'' in your MATLAB path!\n');
            return;
        else
            flag = find(codeDir=='\' | codeDir=='/',1,'last');
            if ~isempty(flag),
                codeDir = codeDir(1 : flag(end)-1);
            end;
            fprintf(1,'\n\nDisplay all STL files in path of stlTools:\n');
            dir([codeDir '/*.stl']);
            if isempty(dir([codeDir '/*.stl'])),
                fprintf(1,'\nNo STL model was found in the path of ''stlTools''!\n');
                return;
            end;
        end;

        fprintf(1,'\nPlease input which type of model nomalization to perform ...\n');
        flag = input('Nomalize along the long axis or all 3 axes? ([]=long axis, other=all 3 axes) ','s');
        norm_sw = isempty(flag);
    end;

    while ~load_stl,
        save_name= input('Name of body model (without stl suffix): ( []=''bee_body'' ) ','s');
        if isempty(save_name),
            save_name = 'bee_body';
        end;
        if ~exist([codeDir '/' save_name '.stl'],'file'),
            fprintf(1,'\nNo file %s found in the path! Please input again!\n', save_name);
            continue;
        end;
        % load and normalize model of body and wing(ascii or binary STL)
        [XYZb,faceb] = stlRead([codeDir '/' save_name '.stl']);
        if isempty(XYZb) || isempty(faceb),
            fprintf(1,'\nSTL file of body found empty! Please input again!\n');
            continue;
        end;
        nptb = size(XYZb,1);
        if norm_sw,
            XYZb =2*(XYZb/(max(XYZb(:,1))-min(XYZb(:,1))))';         % normalized by body length (radius to 1)
        else
            XYZb = 2*(XYZb./(ones(nptb,1)*(max(XYZb,[],1)-min(XYZb, [],1))))';         % normalized like cube (radius to 1)
        end;

        if func_flag,
            disp('Please input name of specific wing model with origin at the root:');
            save_name= input('Name of wing model (without stl suffix): ( []=''bee_wingr'' ) ','s');
            if isempty(save_name),
                save_name = 'bee_wingr';
            end;
        else
            disp('Please input name of specific wing model with origin at the center:');
            save_name= input('Name of wing model (without stl suffix): ( []=''bee_wingc'' ) ','s');
            if isempty(save_name),
                save_name = 'bee_wingc';
            end;
        end;
        if ~exist([codeDir '/' save_name '.stl'],'file'),
            fprintf(1,'\nNo file %s found in the path! Please input again!\n', save_name);
            continue;
        end;
        [XYZw,facew] = stlRead([codeDir '/' save_name '.stl']);
        if isempty(XYZw) || isempty(facew),
            fprintf(1,'\nSTL file of wing found empty! Please input again!\n');
            continue;
        end;
        nptw = size(XYZw,1);
        if norm_sw,
            %     XYZw =2*(XYZw/max(XYZw(:,2)))';        % normalized by wing length from origin to tip (radius to 1)
            XYZw =2*(XYZw/(max(XYZw(:,2))-min(XYZw(:,2))))';        % normalized by wing length (radius to 1)
        else
            % normalized like cube (radius to 1)
            XYZw = 2*(XYZw./(ones(nptw,1)*(max(XYZw,[],1)-min(XYZw, [],1))))';
        end;
        load_stl = 1;
    end;
    % call function
    if func_flag,
        % compute wing root position in body system first
        polygon_projection;
    else
        polygon_center_projection;
    end;

else

    % mesh of unit sphere
    Nb = 10;
    npts = (Nb+1)^2;
    [Xe, Ye, Ze] = sphere(Nb);
    [faces,XYZs] = surf2patch(Xe,Ye,Ze);
    XYZs = XYZs';
    nn = npts+3;
    Nb = nn*(nparts+1);
    % call function
    ellipsoid_projection;
end;

fprintf(1,'\nDone with the reprojection process...\n\nCheck reprojected images in your directory...\n');
%% prepare for reconstruction from reprojection
if strcmp(fgprefix, 'silhouette_'),
    save_name = [imgdir '/old/kinematic_animation.mat'];
    assert(exist(save_name, 'file')==2, 'Old ''kinematic_animation.mat'' not found in ''old'' directory!');
    return;
end;
fprintf(1,'\nYou can reconstruct visual hull and extract kinematics again from reprojection!\n');
flag = input('Prepare to reconstruct visual hull from reprojection or not? ([]=yes, other=no) ','s');
if isempty(flag),
    filepath = [imgdir '/old'];
    if ~exist(filepath,'dir'),
        mkdir(imgdir,'old');
    end;
    disp('Move current visual hull and kinematic data file to directory ''old''...');
    save_name = ['Extraction_result.mat'];
    movefile([imgdir '/' save_name],[filepath '/' save_name]);
    save_name = 'visualhull_environment.mat';
    movefile([imgdir '/' save_name],[filepath '/' save_name]);
    for ii=1:level,
        save_name = ['visual_hull_level' num2str(ii) '.mat'];
        movefile([imgdir '/' save_name],[filepath '/' save_name]);
    end;
    for nl = level : 1+n_sublevel,
        for ii=1:n_unit,
            save_name = ['level' num2str(nl) '_part' sprintf(['%0' ndigit 'd'],ii) '.mat'];
            movefile([imgdir '/' save_name],[filepath '/' save_name]);
        end;
    end;
    for ii=1:n_unit,
        nc = sprintf(['%0' ndigit 'd'],ii);
        save_name = ['3D_points_part' nc '.mat'];
        movefile([imgdir '/' save_name],[filepath '/' save_name]);
        save_name = ['segment_part' nc '.mat'];
        movefile([imgdir '/' save_name],[filepath '/' save_name]);
    end;
    save_name = 'segment_variables.mat';
    movefile([imgdir '/' save_name],[filepath '/' save_name]);
    save_name = ['kinematic_animation.mat'];
    movefile([imgdir '/' save_name],[filepath '/' save_name]);
    save_name = ['wing_axes.jpg'];
    movefile([imgdir '/' save_name],[filepath '/' save_name]);
    save_name = ['wing_axes.fig'];
    movefile([imgdir '/' save_name],[filepath '/' save_name]);
    for ii=0:1,
        save_name = ['kinematics' num2str(ii) '.jpg'];
        movefile([imgdir '/' save_name],[filepath '/' save_name]);
        save_name = ['kinematics' num2str(ii) '.fig'];
        movefile([imgdir '/' save_name],[filepath '/' save_name]);
    end;
    save_name = ['Av_eulerAngle.jpg'];
    movefile([imgdir '/' save_name],[filepath '/' save_name]);
    save_name = ['Av_eulerAngle.fig'];
    movefile([imgdir '/' save_name],[filepath '/' save_name]);
    if calib_mode,
        for kk = ind_active,
            framenb = strnum_frame{kk};
            save_name = ['projection_' imgbase framenb '.jpg'];
            movefile([imgdir '/' save_name],[filepath '/' save_name]);
            save_name = ['3D_projection_' framenb '.jpg'];
            movefile([imgdir '/' save_name],[filepath '/' save_name]);
            save_name = ['3D_projection_' framenb '.fig'];
            movefile([imgdir '/' save_name],[filepath '/' save_name]);
        end;
    else
        for pp = 1:n_cam,
            mkdir(filepath,['Cam_' num2str(pp)]);
        end;
        for kk = ind_active,
            for pp = 1:n_cam,
                if active_imgviews(pp,kk),
                    save_name = ['projection_' imgbase_cell{pp} strnum_frameNcam{pp}{kk} '.jpg'];
                    movefile([imgdir_cell{pp} '/' save_name],[filepath '/Cam_' num2str(pp) '/' save_name]);
                end;
            end;
            framenb = strnum_frameNcam{1}{kk};
            save_name = ['3D_projection_' framenb '.jpg'];
            movefile([imgdir '/' save_name],[filepath '/' save_name]);
            save_name = ['3D_projection_' framenb '.fig'];
            movefile([imgdir '/' save_name],[filepath '/' save_name]);
        end;
    end;

    % update variables for extraction result
    fgprefix = 'silhouette_';
    save_name = [imgdir '/Extraction_result.mat'];
    string_save = ['save ' save_name ' fgprefix calib_mode active_imgviews n_frame n_cam handcc Omcc Tcc'];
    if exist('Omcc_error','var'),
        string_save = [string_save ' Omcc_error Tcc_error'];
    end;
    if calib_mode,
        string_save = [string_save ' imgbase imgfmt strnum_frame nx ny fc cc alpha_c kc'];
        if exist('fc_error','var'),
            string_save = [string_save ' fc_error cc_error alpha_c_error kc_error'];
        end;
    else
        string_save = [string_save ' imgbase_cell imgfmt_cell strnum_frameNcam' ...
            ' imsize fc_mat cc_mat alpha_vec kc_mat'];
        if exist('fc_mat_error','var'),
            string_save = [string_save ' fc_mat_error cc_mat_error alpha_vec_error kc_mat_error'];
        end;
    end;
    eval(string_save);
    fprintf(1,'\nDone...\nYou can reconstruct visualhull now!\n');
end;
