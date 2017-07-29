% 多摄像机参数导入和存储
% Omcc,Tcc,handcc是所有摄像机相对于第一个摄像机的运动变换，所以第一列为cam1到自己的变换

if exist('Camera_list.mat','file')==2,
    fprintf(1,'\nLoading camera parameters from ''Camera_list.mat''...\n');
    load('Camera_list.mat');
    return;
end;
fprintf(1,'\nCamera parameters file ''Camera_list.mat'' not found!\n');
fprintf(1,['\nDo you want to load camera parameters from calibration results ([]=yes=1),\n ' ...
    'or do you want to input them manually (no=0)?\n']);
load_calib = input('Do you want to load camera parameters automatically? ([]=yes=1)');
if isempty(load_calib),
    load_calib = 1;
else
    load_calib = ~~load_calib;
end;

dir('Multi*_Calib_Results*');        % display all calibration results of mirror and multicam
fprintf(1,['\nWhich configuration of calibration result do you want to use?\n' ...
    '(Mode: One camera with split views; Multiple different cameras) ...\n']);
calib_mode = input('calib_mode = ([]=one camera, other=different cameras) ','s');
calib_mode = isempty(calib_mode);

if load_calib,
    fprintf(1,'\nWhich format of calibration result file do you want to load (suffix: m, mat)?\n');
    calib_format = input('calib_format = ([]=m=''m'', other=''mat'') ' ,'s');
    if isempty(calib_format) || strcmp(calib_format,'m'),
        calib_format = 'm';
    else
        calib_format = 'mat';
    end;
    if calib_mode,         % load (one camera with split views) calibration result file
        if strcmp(calib_format,'m'),
            if exist('Multimirror_Calib_Results.m','file')==2,
                Multimirror_Calib_Results;
            else
                fprintf(1,'\nMultimirror_Calib_Results.m not found! Please input camera parameters manually!\n');
                load_calib = 0;
            end;
        else
            if exist('Multimirror_Calib_Results.mat','file')==2,
                load('Multimirror_Calib_Results.mat');
            else
                fprintf(1,'\nMultimirror_Calib_Results.mat not found! Please input camera parameters manually!\n');
                load_calib = 0;
            end;
        end;
    else                        % load (multiple different cameras) calibration result file
        if strcmp(calib_format,'m'),
            if exist('Multicam_Calib_Results.m','file')==2,
                Multicam_Calib_Results;
            else
                fprintf(1,'\nMulticam_Calib_Results.m not found! Please input camera parameters manually!\n');
                load_calib = 0;
            end;
        else
            if exist('Multicam_Calib_Results.mat','file')==2,
                load('Multicam_Calib_Results.mat');
            else
                fprintf(1,'\nMulticam_Calib_Results.mat not found! Please input camera parameters manually!\n');
                load_calib = 0;
            end;
        end;
    end;
end;

if exist('n_cam','var') && n_cam<2,
    fprintf(1,['\nIt takes at least 2 cameras for 3D reconstruction!\n' ...
        'You have to input camera parameters manually!\n']);
    clear n_cam hand_list Omcc Tcc Omcc_error Tcc_error nx ny fc cc alpha_c kc ...
        fc_error cc_error alpha_c_error kc_error imsize fc_mat cc_mat alpha_vec ...
        kc_mat fc_mat_error cc_mat_error alpha_vec_error kc_mat_error;
    load_calib = 0;
end;

% input camera parameters manually
if ~load_calib,
    fprintf(1,'\nManually input intrinsic and extrinsic parameters for all cameras ... \n');
    flag = ~exist('n_cam','var');
    while flag,
        fprintf(1,'\nPlease input number of camera views in your configuration ...\n');
        n_cam = input('Number of cameras = ([]=3) ');
        if isempty(n_cam),
            n_cam = 3;
        end;
        flag = n_cam<2;
        if flag,
            fprintf(1,'\nIt takes at least 2 cameras for 3D reconstruction! Please input again!\n');
        end;
    end;
    % hand_list(i) means axis handness of the i-th camera
    if ~exist('hand_list','var'),
        flag = 1;
    else
        hand_list =  hand_list(:)';
        flag = ~isequal(abs(hand_list),ones(1,n_cam));
        if flag,
            fprintf(1,'\nUnexpected value for handedness vector! Please input again!\n');
        end;
    end;
    while flag,
        fprintf(1,['\nPlease input the axis handedness vector of all %d camera views!\n' ...
            'Note: 1 stands for right-handed, -1 stands for left-handed.\n'],n_cam);
        hand_list = input(['handedness vector = ([] = [' num2str(ones(1,n_cam)) '])']);
        if isempty(hand_list),
            hand_list = ones(1,n_cam);
            flag = 0;
        else
            hand_list =  hand_list(:)';
            flag = ~isequal(abs(hand_list),ones(1,n_cam));
            if flag,
                fprintf(1,'\nUnexpected value for handedness vector! Please input again!\n');
            end;
        end;
    end;
    handcc = hand_list(1)*hand_list;
    
    % Omcc(:,i) is the relative orientation of the i-th camera wrt camera 1 (or reference frame)
    if ~exist('Omcc','var') || ~exist('Tcc','var'),
        flag = 1;
    else
        flag = ~(isequal(size(Omcc),[3,n_cam]) & isequal(size(Tcc),[3,n_cam]));
        if flag,
            fprintf(1,'\nUnexpected dimension for axis angle rotation or translation! Please input again!\n');
        end;
    end;
    if flag,
        fprintf(1,'\nPlease input the axis angle rotation and translation of all %d cameras!\n',n_cam);
        Omcc = zeros(3,n_cam);
        Tcc = Omcc;
        for pp=1:n_cam,
            Omcc(:,pp) = input(['The axis angle orientation vector of the of camera ', num2str(pp) ' = ']);
            Tcc(:,pp) = input(['The translation vector of the of camera ', num2str(pp) ' = ']);
        end;
    end;
    
    if calib_mode,          % one camera with mirrors
        if ~exist('nx','var') || exist('ny','var'),
            nx = input('Width of image = ');
            ny = input('Height of image = ');
        end;
        if ~exist('fc','var'),
            fc = input('Focal length in pixel = ([] = [1;1])');
            if isempty(fc),
                fc = [1;1];
            end;
        end;
        if ~exist('cc','var'),
            cc = input('Principle point in pixel = ([] = ([nx;ny]-1)/2)');
            if isempty(cc),
                cc = ([nx;ny]-1)/2;
            end;
        end;
        if ~exist('alpha_c','var'),
            alpha_c = input('Skew coefficient = ([] = 0)');
            if isempty(alpha_c),
                alpha_c = 0;
            end;
        end;
        if ~exist('kc','var'),
            kc = input('Distortion coefficients = ([] = [0;0;0;0;0])');
            if isempty(kc),
                kc = zeros(5,1);
            end;
        end;
    else                        % multiple different camera
        if ~exist('imsize','var') || ~isequal(size(imsize),[2,n_cam]),
            imsize = ones(2,n_cam);
            for pp = 1:n_cam,
                imsize(:,pp) = input(['Image size [nx; ny] of camera ' num2str(pp) ' = ']);
            end;
        end;
        if ~exist('fc_mat','var') || ~isequal(size(fc_mat),[2,n_cam]),
            fc_mat = ones(2,n_cam);
            for pp = 1:n_cam,
                fc = input(['Focal length of camera ' num2str(pp) ' in pixel = ([] = [1;1])']);
                if isempty(fc),
                    fc = [1;1];
                end;
                fc_mat(:,pp) = fc;
            end;
        end;
        if ~exist('cc_mat','var') || ~isequal(size(cc_mat),[2,n_cam]),
            cc_mat = (imsize-1)/2;
            for pp = 1:n_cam,
                cc = input(['Principle point of camera ' num2str(pp) ' in pixel = ([] = ([nx;ny]-1)/2)']);
                if ~isempty(cc),
                    cc_mat(:,pp) = cc;
                end;
            end;
        end;
        if ~exist('alpha_vec','var'),
            flag = 1;
        else
            alpha_vec = alpha_vec(:)';
            flag = length(alpha_vec)~=n_cam;
            if flag,
                fprintf(1,'\nDimension of alpha_vec do not match number of cameras! Please input again!\n');
            end;
        end;
        while flag,
            alpha_vec = input(['alpha_vec = ([] = [' num2str(zeros(1,n_cam)) '])']);
            if isempty(alpha_vec),
                alpha_vec = zeros(1,n_cam);
                flag = 0;
            else
                alpha_vec = alpha_vec(:)';
                flag = length(alpha_vec)~=n_cam;
                if flag,
                    fprintf(1,'\nDimension of alpha_vec do not match number of cameras! Please input again!\n');
                end;
            end;
        end;
        if ~exist('kc_mat','var') || ~isequal(size(kc_mat),[5,n_cam]),
            kc_mat = zeros(5,n_cam);
            for pp = 1:n_cam,
                kc = input(['Distortion coefficients of camera ' num2str(pp) ' = ([] = [0;0;0;0;0])']);
                if isempty(kc),
                    kc = zeros(5,1);
                end;
                kc_mat(:,pp) = kc;
            end;
        end;
    end;
end;

string_save = 'save Camera_list calib_mode n_cam handcc Omcc Tcc';
if exist('Omcc_error','var'),
    string_save = [string_save ' Omcc_error Tcc_error'];
end;
if calib_mode,
    string_save = [string_save ' nx ny fc cc alpha_c kc'];
    if exist('fc_error','var'),
        string_save = [string_save ' fc_error cc_error alpha_c_error kc_error'];
    end;
else
    string_save = [string_save ' imsize fc_mat cc_mat alpha_vec kc_mat'];
    if exist('fc_mat_error','var'),
        string_save = [string_save ' fc_mat_error cc_mat_error alpha_vec_error kc_mat_error'];
    end;
end;
eval(string_save);
fprintf(1,'done...\nCamera parameters is stored in ''Camera_list.mat''.\n');