% convert the classical Bouguet format to my format (Multiple different cameras) or vice versa
% you can convert multiple calib_data*.mat from different camera views of same pattern to
% multicam_calib_data.mat or vice versa

fprintf(1,['\nYou can convert the classical Bouguet format to multiple-different-cameras\n' ...
    'format (calib_data to multicam_calib_data) or choose the other way around ...\n']);
conv_flag = input('Please choose alternative ([]=Bouguet to multicam, other=vice versa) ', 's');
if isempty(conv_flag),
    fprintf(1,'\nConvert files ''calib_data*.mat'' to ''multicam_calib_data.mat''.\n');
    dir('calib_data*.mat');
    datasets = dir('calib_data*.mat');
    Ndata = length(datasets);  % number of calib_data*.mat
    if Ndata==0,
        fprintf(1,'No file ''calib_data*.mat'' in current directory!\n');
        return;
    else
        % load all variables of the 1st dataset seperately into workspace
        save_name = datasets(1).name;
        load(save_name);
        n_cam = Ndata;
        if n_cam>1,
            fprintf(1,'Make sure all ''calib_data*.mat'' files contain corresponding points of multiple cameras!\n');
            flag = input('Do these files constitute different camera views of same pattern? ([]=no, other=yes)');
            if isempty(flag),
                disp('Please remove extra data to meet the requirement!')
                return;
            end;
            % check data
            for pp = 2:n_cam,
                % load variables of dataset pp into cell datapp
                datapp = load(datasets(pp).name);
                if n_ima~=datapp.n_ima,
                    fprintf(1,'\nERROR: Number of images in data #%d do not match with data #1!\n',pp);
                    return;
                end;
            end;
        end;
        
        % initialization
        % input handedness vector for every camera
        flag = 1;
        while flag,
            fprintf(1,['Please input the axis handedness vector of all ' num2str(n_cam) ' cameras!\n' ...
                'Note: 1 stands for right-handed, -1 stands for left-handed.\n']);
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
        active_imgviews = zeros(n_cam,n_ima);
        n_view = n_ima * n_cam;
        X_cell =  cell(1,n_view);
        x_cell =  X_cell;
        dXY_mat = NaN(2,n_view);
        n_sq_mat =  dXY_mat;
        imbase = cell(1,n_cam);
        imformat = imbase;
        imstrnum = imbase;
        imsize = zeros(2,n_cam);
        win_size = imsize;       
        % load data for all views
        for pp = 1:n_cam,
            % load all variables of dataset pp seperately into workspace
            save_name = datasets(pp).name;
            fprintf(1,'Load file ''%s'' ...\n',save_name);
            load(save_name);
            strnum_cell = cell(1,n_ima);
            hand = hand_list(pp);
            for kk = 1:n_ima,
                % load string_num
                if ~type_numbering,
                    number_ext =  num2str(image_numbers(kk));
                else
                    number_ext = sprintf(['%0' num2str(N_slots) 'd'],image_numbers(kk));
                end;
                strnum_cell{kk} = number_ext;
                X_kk = eval(['X_' num2str(kk)]);
                x_kk = eval(['x_' num2str(kk)]);
                dX = eval(['dX_' num2str(kk)]);
                dY = eval(['dY_' num2str(kk)]);
                n_sq_x = eval(['n_sq_x_' num2str(kk)]);
                n_sq_y = eval(['n_sq_y_' num2str(kk)]);
                % transform right hand frame to left hand (symmetric operation)
                if hand~=1,
                    % initialize array column index
                    % four corners in counterclockwise direction (o, x, xy, y)
                    idx = 1:(n_sq_x+1)*(n_sq_y+1);
                    idx = reshape(idx, n_sq_x+1, n_sq_y+1);
                    % normal grid, four corners in counterclockwise direction (y,xy,x,o)
                    idx = fliplr(idx);
                    % rotate array 90 degree clockwise to get left hand frame (swap x and y axis)
                    idx = rot90(idx,-1);
                    idx = idx(:)';                       % column index to swap x and y axis
                    X_kk = X_kk([2 1 3],idx);   % swap X and Y to finish transformation
                    x_kk = x_kk(:,idx);         % pixel points for split views only change index
                    % swap dX and dY
                    temp = dX;
                    dX = dY;
                    dY = temp;
                    % swap n_sq_x and n_sq_y
                    temp = n_sq_x;
                    n_sq_x = n_sq_y;
                    n_sq_y = temp;
                end;
                kth = (kk-1)*n_cam+pp;
                X_cell{kth} = X_kk;
                x_cell{kth} = x_kk;
                dXY_mat(1,kth) = dX;
                dXY_mat(2,kth) = dY;
                n_sq_mat(1,kth) = n_sq_x;
                n_sq_mat(2,kth) = n_sq_y;
            end;
            active_imgviews(pp,:) = active_images;
            imsize(1,pp) = nx;
            imsize(2,pp) = ny;
            win_size(1,pp) = wintx;
            win_size(2,pp) = winty;
            imbase{pp} = calib_name;
            imformat{pp} = format_image;
            imstrnum{pp} = strnum_cell;
        end;
        active_images = any(active_imgviews,1);
        ind_active = find(active_images);
        check_extracted_mirror; % Fix nan variables to []
        
        % check save_name
        save_name = 'multicam_calib_data';
        if exist([save_name '.mat'],'file')==2,
            disp('WARNING: File ''multicam_calib_data.mat'' already exists!');
            if exist('copyfile','builtin')==5,
                cont = -1;
                flag = 1;
                while flag,
                    cont = cont + 1;
                    save_name = ['old_multicam_calib_data' num2str(cont)];
                    flag = (exist([save_name '.mat'],'file')==2);
                end;
                copyfile('multicam_calib_data.mat',[save_name '.mat']);
                disp(['Copying the current ''multicam_calib_data.mat'' file to ''' save_name '.mat''']);
                cont_save = 1;
            else
                disp('The file multicam_calib_data.mat is about to be changed.');
                cont_save = input('Do you want to continue? ([]=no,other=yes) ','s');
                cont_save = ~isempty(cont_save);
            end;
        else
            cont_save = 1;
        end;
        if ~cont_save,
            fprintf(1,'\nNew format calibration data are not saved!\n');
            return;
        end;
        save_name = 'multicam_calib_data';
        string_save = ['save ' save_name ' active_images ind_active wintx winty n_ima map ' ...
            'dX_default dY_default dX dY X_cell x_cell n_sq_mat dXY_mat win_size n_cam ' ...
            'active_imgviews hand_list imsize imbase imformat imstrnum'];
        eval(string_save);
        disp('Done! New format of data are saved in ''multicam_calib_data.mat''.');
    end;
    
else
    %% convert my format (one camera with split views) to the classical Bouguet format
    fprintf(1,'\nConvert file ''multicam_calib_data.mat'' to separate ''calib_data*.mat'' files.\n');
    dir('multicam_calib_data*.mat');
    save_name = 'multicam_calib_data';
    if exist([save_name '.mat'],'file')~=2,
        fprintf(1,'File ''multicam_calib_data.mat'' do not exist!\n');
        return;
    end;
    datapp = load(save_name);
    % number of cameras determine number of calib_data*.mat
    n_cam = datapp.n_cam;
    % number of valid images
    Nima_valid = datapp.n_ima;
    map = datapp.map;
    dX_default = datapp.dX_default;
    dY_default = datapp.dY_default;
    
    % my format do not need image number to be continuous
    % while Bouguet format do. See check_directory.m for detail.
    % load data from camera 1
    calib_name = datapp.imbase{1};
    format_image = datapp.imformat{1};
    strnum_cell = datapp.imstrnum{1};
    indices =  zeros(1,Nima_valid);
    string_length = indices;
    for kk = 1:Nima_valid,
        string_num = strnum_cell{kk};
        string_length(kk) = size(string_num,2);
        indices(kk) = str2double(string_num);
    end;
    first_num = indices(1);
    n_ima = indices(end) - first_num + 1;
    if min(string_length) == max(string_length),
        N_slots = min(string_length);
        type_numbering = 1;
    else
        N_slots = 1;
        type_numbering = 0;
    end;
    image_numbers = first_num : n_ima-1+first_num;

    nx = datapp.imsize(1,1);
    ny = datapp.imsize(2,1);
    Wcal = nx;
    Hcal = ny;
    wintx = datapp.win_size(1,1);
    winty = datapp.win_size(2,1);
    active_view = datapp.active_imgviews(1,:);
    active_images = zeros(1,n_ima);
    hand = datapp.hand_list(1);
    for kk = 1:Nima_valid,
        kki = find(image_numbers==indices(kk));
        if active_view(kk),
            active_images(kki) = 1;
            kth = (kk-1)*n_cam+1;
            X_kk = datapp.X_cell{kth};
            [m, n] = size(X_kk);
            if m==2,
                X_kk = [X_kk; zeros(1,n)];
            end;
            x_kk = datapp.x_cell{kth};
            dX = datapp.dXY_mat(1,kth);
            dY = datapp.dXY_mat(2,kth);
            n_sq_x = datapp.n_sq_mat(1,kth);
            n_sq_y = datapp.n_sq_mat(2,kth);
            % transform left hand frame back to right hand (symmetric operation)
            if hand~=1,
                % swap x and y axis
                idx = 1:(n_sq_x+1)*(n_sq_y+1);
                idx = reshape(idx, n_sq_x+1, n_sq_y+1);
                idx = fliplr(idx);
                idx = rot90(idx,-1);
                idx = idx(:)';
                X_kk = X_kk([2 1 3],idx);   % swap X and Y to finish transformation
                x_kk = x_kk(:,idx);       % pixel points for split views only change index
                % swap dX and dY
                temp = dX;
                dX = dY;
                dY = temp;
                % swap n_sq_x and n_sq_y
                temp = n_sq_x;
                n_sq_x = n_sq_y;
                n_sq_y = temp;
            end;
            eval(['X_' num2str(kki) ' = X_kk;']);
            eval(['x_' num2str(kki) ' = x_kk;']);
            eval(['dX_' num2str(kki) ' = dX;']);
            eval(['dY_' num2str(kki) ' = dY;']);
            eval(['n_sq_x_' num2str(kki) ' = n_sq_x;']);
            eval(['n_sq_y_' num2str(kki) ' = n_sq_y;']);
            eval(['wintx_' num2str(kki) ' = wintx;']);
            eval(['winty_' num2str(kki) ' = winty;']);
        end;
    end;
    ind_active = find(active_images);
    
    % Fix potential non-existing variables:
    for kk = 1:n_ima,
        if ~exist(['x_' num2str(kk)], 'var'),
            eval(['X_' num2str(kk) ' = NaN*ones(3,1);']);
            eval(['x_' num2str(kk) ' = NaN*ones(2,1);']);
            eval(['dX_' num2str(kk) ' = NaN;']);
            eval(['dY_' num2str(kk) ' = NaN;']);
            eval(['n_sq_x_' num2str(kk) ' = NaN;']);
            eval(['n_sq_y_' num2str(kk) ' = NaN;']);
        end;
        if ~exist(['wintx_' num2str(kk)],'var') || ~exist(['winty_' num2str(kk)],'var'),
            eval(['wintx_' num2str(kk) ' = NaN;']);
            eval(['winty_' num2str(kk) ' = NaN;']);
        end;
    end;
    
    % save new data: only one camera
    if n_cam==1,
        % check save_name
        save_name = 'calib_data';
        if exist([save_name '.mat'],'file')==2,
            disp('WARNING: File ''calib_data.mat'' already exists!');
            if exist('copyfile','builtin')==5,
                cont = -1;
                flag = 1;
                while flag,
                    cont = cont + 1;
                    save_name = ['old_calib_data' num2str(cont)];
                    flag = (exist([save_name '.mat'],'file')==2);
                end;
                copyfile('calib_data.mat',[save_name '.mat']);
                disp(['Copying the current ''calib_data.mat'' file to ''' save_name '.mat''']);
                cont_save = 1;
            else
                disp('The file ''calib_data.mat'' is about to be changed.');
                cont_save = input('Do you want to continue? ([]=no,other=yes) ','s');
                cont_save = ~isempty(cont_save);
            end;
        else
            cont_save = 1;
        end;
        if ~cont_save,
            fprintf(1,'\nNew format calibration data are not saved!\n');
            return;
        end;
        save_name = 'calib_data';
        string_save = ['save ' save_name ' active_images ind_active wintx winty n_ima type_numbering N_slots ' ...
            'first_num image_numbers format_image calib_name Hcal Wcal nx ny map dX_default dY_default dX dY'];
        for kk = 1:n_ima,
            string_save = [string_save ' X_' num2str(kk) ' x_' num2str(kk) ' n_sq_x_' num2str(kk) ' n_sq_y_' num2str(kk) ...
                ' wintx_' num2str(kk) ' winty_' num2str(kk) ' dX_' num2str(kk) ' dY_' num2str(kk)];
        end;
        eval(string_save);
        disp('Done! New format of data are saved in ''calib_data.mat''.');
        
    else
        % multiple cameras
        fprintf(1,['\nNew format of data will be stored in ''calib_data(1~%d).mat''.\n' ...
            'File with same name will be erased!\n'],n_cam);
        flag = input('Are you sure to save these files? ([]=no, other=yes) ','s');
        if isempty(flag),
            fprintf(1,'\nNew format calibration data are not saved!\n');
            return;
        end;
        % save camera 1
        save_name = 'calib_data1';
        string_save = ['save ' save_name ' active_images ind_active wintx winty n_ima type_numbering N_slots ' ...
            'first_num image_numbers format_image calib_name Hcal Wcal nx ny map dX_default dY_default dX dY'];
        for kk = 1:n_ima,
            string_save = [string_save ' X_' num2str(kk) ' x_' num2str(kk) ' n_sq_x_' num2str(kk) ' n_sq_y_' num2str(kk) ...
                ' wintx_' num2str(kk) ' winty_' num2str(kk) ' dX_' num2str(kk) ' dY_' num2str(kk)];
        end;
        eval(string_save);
        
        % save camera 2 and above
        for pp = 2:n_cam,
            calib_name = datapp.imbase{pp};
            format_image = datapp.imformat{pp};
            strnum_cell = datapp.imstrnum{pp};
            indices =  zeros(1,Nima_valid);
            string_length = indices;
            for kk = 1:Nima_valid,
                string_num = strnum_cell{kk};
                string_length(kk) = size(string_num,2);
                indices(kk) = str2double(string_num);
            end;
            first_num = indices(1);
            % n_ima should be same as camera 1, but N_slots may be different
            n_ima = indices(end) - first_num + 1;
            if min(string_length) == max(string_length),
                N_slots = min(string_length);    % digits of image number
                type_numbering = 1;
            else
                N_slots = 1;
                type_numbering = 0;
            end;
            image_numbers = first_num : n_ima-1+first_num;
            
            nx = datapp.imsize(1,pp);
            ny = datapp.imsize(2,pp);
            Wcal = nx;
            Hcal = ny;
            wintx = datapp.win_size(1,pp);
            winty = datapp.win_size(2,pp);
            active_view = datapp.active_imgviews(pp,:);
            active_images = zeros(1,n_ima);
            hand = datapp.hand_list(pp);
            for kk = 1:Nima_valid,
                kki = find(image_numbers==indices(kk));
                if active_view(kk),
                    active_images(kki) = 1;
                    kth = (kk-1)*n_cam+pp;
                    X_kk = datapp.X_cell{kth};
                    [m, n] = size(X_kk);
                    if m==2,
                        X_kk = [X_kk; zeros(1,n)];
                    end;
                    x_kk = datapp.x_cell{kth};
                    dX = datapp.dXY_mat(1,kth);
                    dY = datapp.dXY_mat(2,kth);
                    n_sq_x = datapp.n_sq_mat(1,kth);
                    n_sq_y = datapp.n_sq_mat(2,kth);
                    % transform left hand frame back to right hand (symmetric operation)
                    if hand~=1,
                        % swap x and y axis
                        idx = 1:(n_sq_x+1)*(n_sq_y+1);
                        idx = reshape(idx, n_sq_x+1, n_sq_y+1);
                        idx = fliplr(idx);
                        idx = rot90(idx,-1);
                        idx = idx(:)';
                        X_kk = X_kk([2 1 3],idx);   % swap X and Y to finish transformation
                        x_kk = x_kk(:,idx);       % pixel points for split views only change index
                        % swap dX and dY
                        temp = dX;
                        dX = dY;
                        dY = temp;
                        % swap n_sq_x and n_sq_y
                        temp = n_sq_x;
                        n_sq_x = n_sq_y;
                        n_sq_y = temp;
                    end;
                    eval(['X_' num2str(kki) ' = X_kk;']);
                    eval(['x_' num2str(kki) ' = x_kk;']);
                    eval(['dX_' num2str(kki) ' = dX;']);
                    eval(['dY_' num2str(kki) ' = dY;']);
                    eval(['n_sq_x_' num2str(kki) ' = n_sq_x;']);
                    eval(['n_sq_y_' num2str(kki) ' = n_sq_y;']);
                    eval(['wintx_' num2str(kki) ' = wintx;']);
                    eval(['winty_' num2str(kki) ' = winty;']);
                end;
            end;
            ind_active = find(active_images);
            
            % Fix potential non-existing variables:
            for kk = 1:n_ima,
                if ~exist(['x_' num2str(kk)], 'var'),
                    eval(['X_' num2str(kk) ' = NaN*ones(3,1);']);
                    eval(['x_' num2str(kk) ' = NaN*ones(2,1);']);
                    eval(['dX_' num2str(kk) ' = NaN;']);
                    eval(['dY_' num2str(kk) ' = NaN;']);
                    eval(['n_sq_x_' num2str(kk) ' = NaN;']);
                    eval(['n_sq_y_' num2str(kk) ' = NaN;']);
                end;
                if ~exist(['wintx_' num2str(kk)],'var') || ~exist(['winty_' num2str(kk)],'var'),
                    eval(['wintx_' num2str(kk) ' = NaN;']);
                    eval(['winty_' num2str(kk) ' = NaN;']);
                end;
            end;
            
            save_name = ['calib_data' num2str(pp)];
            string_save = ['save ' save_name ' active_images ind_active wintx winty n_ima type_numbering N_slots ' ...
                'first_num image_numbers format_image calib_name Hcal Wcal nx ny map dX_default dY_default dX dY'];
            for kk = 1:n_ima,
                string_save = [string_save ' X_' num2str(kk) ' x_' num2str(kk) ' n_sq_x_' num2str(kk) ' n_sq_y_' num2str(kk) ...
                    ' wintx_' num2str(kk) ' winty_' num2str(kk) ' dX_' num2str(kk) ' dY_' num2str(kk)];
            end;
            eval(string_save);
        end;
        fprintf('\nDone! New format of data are saved in ''calib_data(1~%d).mat''.\n',n_cam);
    end;
end;