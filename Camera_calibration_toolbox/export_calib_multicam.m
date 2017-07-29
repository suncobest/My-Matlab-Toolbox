%% Export calibration data (corners + 3D coordinates) to
%% text files (in Willson-Heikkila's format or Zhang's format)

%% Thanks to Michael Goesele (from the Max-Planck-Institut) for the original suggestion
%% of adding this export function to the toolbox.

if ~exist('x_cell', 'var'),
    fprintf(1,'ERROR: No calibration data to export!\n');
    return;
end;

if n_ima == 0,
    fprintf(1,'ERROR: No calibration data to export!\n');
    return;
end;
active_images = any(active_imgviews,1);
ind_active = find(active_images);
for kk =  ind_active,
    for pp = 1:n_cam,
        kth = (kk-1)*n_cam+pp;
        x = x_cell{kth};
        if isempty(x) || isnan(x(1)),   % 若x为[]或NaN，则此视角无效；
            active_imgviews(pp,kk) = 0;
            x_cell{kth} = [];
            X_cell{kth} = [];
            dXY_mat(:, kth) = NaN(2,1);
            n_sq_mat(:, kth) = NaN(2,1);
        end;
    end;
    if all(active_imgviews(:,kk)==0),
        fprintf(1,['WARNING:  No camera have grid corners on image ' ...
            num2str(kk) '- This image is now set inactive!\n']);
    end;
end;
active_images = any(active_imgviews,1);
ind_active = find(active_images);

fprintf(1,'Exports calibration data to Willson-Heikkila or Zhang formats:\n');
fprintf(1,'Two possible formats of export: 0=[]=Willson and Heikkila, 1=Zhang\n');
flag = input('Format of export (enter 0=[] or 1): ');
if isempty(flag)
    flag = 0;
else
    flag = ~~flag;
end;

if ~exist('n_view','var'),
    n_view = n_ima * n_cam;
end;
ndigit = num2str(floor(log10(n_view))+1);
ind_active_views = find(active_imgviews(:)');

if flag,
    fprintf(1,'\nExport of calibration data to text files (Zhang''s format)\n');
    model_base = input('File basename for the 3D world coordinates: ( []=X3d_ )','s');
    if isempty(model_base),
        model_base = 'X3d_';
    end;
    pixel_base = input('File basename for the 2D image coordinates: ( []=x2d_ )','s');
    if isempty(pixel_base),
        pixel_base = 'x2d_';
    end;
    for kth = ind_active_views,
        X_kk = X_cell{kth};
        x_kk = x_cell{kth};

        n_sq_x = n_sq_mat(1,kth);
        n_sq_y = n_sq_mat(2,kth);

        % Matrix start from y to xy along x direction: (x,y)
        X = reshape(X_kk(1,:),n_sq_x+1,n_sq_y+1)';
        Y = reshape(X_kk(2,:),n_sq_x+1,n_sq_y+1)';
        XY = reshape([X;Y],n_sq_y+1,2*(n_sq_x+1));

        x = reshape(x_kk(1,:),n_sq_x+1,n_sq_y+1)';
        y = reshape(x_kk(2,:),n_sq_x+1,n_sq_y+1)';
        xy = reshape([x;y],n_sq_y+1,2*(n_sq_x+1));

        kk = floor((kth-1)/n_cam)+1;
        pp = kth-(kk-1)*n_cam;
        extnum = sprintf(['%0' ndigit 'd'],kth);
        fprintf(1, ['Exporting calibration data of (camera ' num2str(pp) ', image ' num2str(kk) ')'...
            '\nto files ' model_base extnum '.txt and ' pixel_base extnum '.txt ...']);
        save([model_base extnum '.txt'], 'XY', '-ASCII');
        save([pixel_base extnum '.txt'], 'xy', '-ASCII');
    end;
else
    fprintf(1,'\nExport of calibration data to text files (Willson and Heikkila''s format)\n');
    outputfile = input('File basename: ( []=Xx_ )','s');
    if isempty(outputfile),
        outputfile = 'Xx_';
    end;
    for kth = ind_active_views,
        X_kk = X_cell{kth};
        x_kk = x_cell{kth};

        % Vector form: [X_kk',x_kk']
        Xx = [X_kk ; x_kk]';

        kk = floor((kth-1)/n_cam)+1;
        pp = kth-(kk-1)*n_cam;
        extnum = sprintf(['%0' ndigit 'd'],kth);
        file_name = [outputfile extnum '.txt'];
        fprintf(1, ['Exporting calibration data  (3D world + 2D image coordinates) of\n'  ...
            '(camera '  num2str(pp) ', image ' num2str(kk) ') to file ' file_name ' ...']);
        save(file_name, 'Xx', '-ASCII');
    end;
end;

fprintf(1,'Done!\n');
