% convert the classical Bouguet format to my format (one camera with split views) or vice versa
% you can convert multiple calib_data*.mat from different views of one camera to
% multimirror_calib_data.mat or vice versa

fprintf(1,['\nYou can convert the classical Bouguet format to one-camera-with-split-views\n' ...
    'format (calib_data to multimirror_calib_data) or choose the other way around ...\n']);
conv_flag = input('Please choose alternative ([]=Bouguet to split views, other=vice versa) ', 's');
if isempty(conv_flag),
    fprintf(1,'\nConvert files ''calib_data*.mat'' to ''multimirror_calib_data.mat''.\n');
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
            fprintf(1,['Make sure all ''calib_data*.mat'' files contain corresponding\n' ...
                'points from different views of one camera!\n']);
            flag = input('All ''calib_data*.mat'' files constitute split views of one camera? ([]=no, other=yes)');
            if isempty(flag),
                disp('Please remove extra data to meet the requirement!')
                return;
            end;
            
            % check data
            for pp = 2:n_cam,
                save_name = datasets(pp).name;
                % load variables of dataset pp into cell datapp
                datapp = load(save_name);
                if n_ima~=datapp.n_ima,
                    fprintf(1,'\nERROR: Number of images in data #%d do not match with data #1!\n',pp);
                    return;
                end;
                if nx~=datapp.nx || ny~=datapp.ny,
                    fprintf(1,'\nERROR: Image sizes in data #%d do not match with data #1!\n',pp);
                    return;
                end;
                if ~strcmp(calib_name, datapp.calib_name),
                    fprintf(1,'\nERROR: Image basename in data #%d do not match with data #1!\n',pp);
                    return;
                end;
                if ~strcmp(format_image, datapp.format_image),
                    fprintf(1,'\nERROR: Image suffix in data #%d do not match with data #1!\n',pp);
                    return;
                end;
                if type_numbering~=datapp.type_numbering,
                    fprintf(1,'\nERROR: Variable type_numbering in data #%d do not match with data #1!\n',pp);
                    return;
                end;
                if N_slots~=datapp.N_slots,
                    fprintf(1,'\nERROR: Variable N_slots in data #%d do not match with data #1!\n',pp);
                    return;
                end;
                if ~isequal(image_numbers, datapp.image_numbers),
                    fprintf(1,'\nERROR: Variable image_numbers in data #%d do not match with data #1!\n',pp);
                    return;
                end;
                if ~isequal(first_num, datapp.first_num),
                    fprintf(1,'\nERROR: Variable first_num in data #%d do not match with data #1!\n',pp);
                    return;
                end;
            end;
        end;
        
        % initialization
        % input handedness vector for every camera
        flag = 1;
        while flag,
            fprintf(1,['Please input the axis handedness vector of all ' num2str(n_cam) ' views!\n' ...
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
        win_size = zeros(2,n_cam);
        n_view = n_ima * n_cam;
        X_cell =  cell(1,n_view);
        x_cell =  X_cell;
        dXY_mat = NaN(2,n_view);
        n_sq_mat =  dXY_mat;
        strnum_cell = cell(1,n_ima);
        % load string_num
        for kk = 1:n_ima,
            if ~type_numbering,
                number_ext =  num2str(image_numbers(kk));
            else
                number_ext = sprintf(['%0' num2str(N_slots) 'd'],image_numbers(kk));
            end;
            strnum_cell{kk} = number_ext;
        end;
        
        % load data for all views
        for pp = 1:n_cam,
            % load all variables of dataset pp seperately into workspace
            save_name = datasets(pp).name;
            fprintf(1,'Load file ''%s'' ...\n',save_name);
            load(save_name);
            active_imgviews(pp,:) = active_images;
            win_size(1,pp) = wintx;
            win_size(2,pp) = winty;
            hand = hand_list(pp);
            for kk = 1:n_ima,
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
                    % rotate array 90 degree clockwise to get left hand frame in mirror (swap x and y axis)
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
        end;
        active_images = any(active_imgviews,1);
        ind_active = find(active_images);
        check_extracted_mirror; % Fix nan variables to []
        
        % check save_name
        save_name = 'multimirror_calib_data';
        if exist([save_name '.mat'],'file')==2,
            disp('WARNING: File ''multimirror_calib_data.mat'' already exists!');
            if exist('copyfile','builtin')==5,
                cont = -1;
                flag = 1;
                while flag,
                    cont = cont + 1;
                    save_name = ['old_multimirror_calib_data' num2str(cont)];
                    flag = (exist([save_name '.mat'],'file')==2);
                end;
                copyfile('multimirror_calib_data.mat',[save_name '.mat']);
                disp(['Copying the current ''multimirror_calib_data.mat'' file to ''' save_name '.mat''']);
                cont_save = 1;
            else
                disp('The file ''multimirror_calib_data.mat'' is about to be changed.');
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
        save_name = 'multimirror_calib_data';
        string_save = ['save ' save_name ' active_images ind_active wintx winty n_ima format_image ' ...
            'calib_name nx ny map dX_default dY_default dX dY X_cell x_cell n_sq_mat dXY_mat ' ...
            'win_size n_cam active_imgviews hand_list strnum_cell'];
        eval(string_save);
        disp('Done! New format of data are saved in ''multimirror_calib_data.mat''.');
    end;
    
else
    %% convert my format (one camera with split views) to the classical Bouguet format
    fprintf(1,'\nConvert file ''multimirror_calib_data.mat'' to separate ''calib_data*.mat'' files.\n');
    dir('multimirror_calib_data*.mat');
    save_name = 'multimirror_calib_data';
    if exist([save_name '.mat'],'file')~=2,
        fprintf(1,'File ''multimirror_calib_data.mat'' do not exist!\n');
        return;
    end;
    datapp = load(save_name);
    % load image name
    calib_name = datapp.calib_name;
    format_image = datapp.format_image;
    % my format do not need image number to be continuous
    % while Bouguet format do. See check_directory.m for detail.
    Nima_valid = datapp.n_ima;
    indices =  zeros(1,Nima_valid);
    string_length = indices;
    for kk = 1:Nima_valid,
        string_num = datapp.strnum_cell{kk};
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
    
    % number of views determine number of calib_data*.mat
    n_cam = datapp.n_cam;
    % load data from view 1
    nx = datapp.nx;
    ny = datapp.ny;
    Wcal = nx;
    Hcal = ny;
    map = datapp.map;
    dX_default = datapp.dX_default;
    dY_default = datapp.dY_default;
    wintx = datapp.win_size(1,1);
    winty = datapp.win_size(2,1);
    active_view = datapp.active_imgviews(1,:);
    ind_images = zeros(1,Nima_valid);
    active_images = zeros(1,n_ima);
    hand = datapp.hand_list(1);
    for kk = 1:Nima_valid,
        kki = find(image_numbers==indices(kk));
        ind_images(kk) = kki;
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
    
    % save new data: only one view
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
        % multiple views
        fprintf(1,['\nNew format of data will be stored in ''calib_data(1~%d).mat''.\n' ...
            'File with same name will be erased!\n'],n_cam);
        flag = input('Are you sure to save these files? ([]=no, other=yes) ','s');
        if isempty(flag),
            fprintf(1,'\nNew format calibration data are not saved!\n');
            return;
        end;
        % save view1
        save_name = 'calib_data1';
        string_save = ['save ' save_name ' active_images ind_active wintx winty n_ima type_numbering N_slots ' ...
            'first_num image_numbers format_image calib_name Hcal Wcal nx ny map dX_default dY_default dX dY'];
        for kk = 1:n_ima,
            string_save = [string_save ' X_' num2str(kk) ' x_' num2str(kk) ' n_sq_x_' num2str(kk) ' n_sq_y_' num2str(kk) ...
                ' wintx_' num2str(kk) ' winty_' num2str(kk) ' dX_' num2str(kk) ' dY_' num2str(kk)];
        end;
        eval(string_save);
        
        % save view2 and above
        for pp = 2:n_cam,
            wintx = datapp.win_size(1,pp);
            winty = datapp.win_size(2,pp);
            active_view = datapp.active_imgviews(pp,:);
            active_images = zeros(1,n_ima);
            hand = datapp.hand_list(pp);
            for kk = 1:Nima_valid,
                if active_view(kk),
                    kki = ind_images(kk);
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
