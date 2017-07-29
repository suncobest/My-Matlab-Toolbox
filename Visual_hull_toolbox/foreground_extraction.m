%%% This script will extract foreground and store binary images
% rectify images before foreground extraction
fprintf(1,['\nFor periodic motion, it is recommended that \n'...
    'the selected image frames contain exactly whole periods!\n']);

if ~exist('calib_mode','var') || ~exist('n_cam','var') ||~exist('Tcc','var'),
    fprintf(1,['\nPlease locate the directory of camera parameter file ''Camera_list.mat'' ...\n']);
    camfolder = uigetdir;
    save_name = [camfolder '/Camera_list.mat'];
    if exist(save_name,'file')~=2,
        fprintf(1,'\nFile ''Camera_list.mat'' not found!\n');
        return;
    end;
    load(save_name);
end;

fprintf(1,['\nAs nonlinear transformation, distortion can dramatically alter the image geometry!\n' ...
    'So image sequence need to be rectified if there is distortion before any further process ...\n']);
flag = input('Image sequences need rectifying or not? ([]=no, other=yes) ','s');
if ~isempty(flag),
    if calib_mode,
        fprintf(1,'\nLocate image frames for one-camera configuation ...\n');
        if norm(kc)==0,
            fprintf(1,'\nNOTE: The camera seems to have no distortion!\n');
        else
            fprintf(1,'\nStart to rectify image sequence ...\n');
            rectify_image_sequence;
            imgbase = ['rect_' imgbase];   % use rectified image
        end;

    else

        fprintf(1,'\nLocate image frames for multiple-camera configuation ...\n');
        if ~exist('imgbase_cell','var'),
            % file variable (input and output)
            imgdir_cell = cell(1,n_cam);
            imgbase_cell = imgdir_cell;
            imgfmt_cell = imgdir_cell;
        end;
        for pp = 1:n_cam,
            % input
            imgdir = imgdir_cell{pp};
            imgbase = imgbase_cell{pp};
            imgfmt = imgfmt_cell{pp};
            % camera parameters
            fc = fc_mat(:,pp);
            cc = cc_mat(:,pp);
            kc = kc_mat(:,pp);
            alpha_c = alpha_vec(pp);

            %  undistort images
            if norm(kc)==0,
                fprintf(1,'\nNOTE: Camera %d seems to have no distortion!\n',pp);
            else
                fprintf(1,'\nStart to rectify image sequence of camera %d...\n',pp);
                rectify_image_sequence;
                imgbase = ['rect_' imgbase];   % use rectified image
            end;

            % output
            imgdir_cell{pp} = imgdir;
            imgbase_cell{pp} = imgbase;
            imgfmt_cell{pp} = imgfmt;
        end;
    end;
end;

fprintf(1,['\nImage sequence are now assumed to be distortion free,\n' ...
    'background images should also be distortion free ...\n']);
flag = input('Have you prepared background images yet? ([]=yes, other=no) ','s');
if ~isempty(flag),
    fprintf(1,'\nYou may use Photoshop to process background images and save with special prefix ...\n');
    return;
end;

if calib_mode,
    fprintf(1,'\nLocate image frames for one-camera configuation ...\n');
    active_imgviews = false(n_cam, n_frame);
    ind_active = [];
    extract_foreground;   % output: imgdir imgbase imgfmt n_frame, strnum_frame
    if isempty(imgdir), return; end;
    if norm(kc)~=0 && ~strcmp(imgbase(1:5), 'rect_'),
        fprintf(1,'\nERROR: You should rectify image sequence first!\n');
        return;
    end;
    active_imgviews(:,ind_active) = 1;
else
    fprintf(1,'\nLocate image frames for multiple-camera configuation ...\n');
    if ~exist('imgbase_cell','var'),
        % file variable (input and output)
        imgdir_cell = cell(1,n_cam);
        imgbase_cell = imgdir_cell;
        imgfmt_cell = imgdir_cell;
    end;
    strnum_frameNcam = cell(1, n_cam);
    active_imgviews = cell(n_cam, 1);
    if ~exist('nframe_slot','var'),
        nframe_slot = zeros(1,n_cam);
        diff_imnum = cell(1,n_cam);
    end;
    for pp = 1:n_cam,
        fprintf(1,'\nFor camera %d ...\n',pp);
        % input
        imgdir = imgdir_cell{pp};
        imgbase = imgbase_cell{pp};
        imgfmt = imgfmt_cell{pp};
        % set background image name and format to empty to allow
        % different names for different cameras
        bgbase = [];
        bgfmt = [];
        ind_active = [];
        % extracting script
        extract_foreground;
        if isempty(imgdir), return; end;
        if norm(kc_mat(:,pp))~=0 && ~strcmp(imgbase(1:5), 'rect_'),
            fprintf(1,'\nERROR: You should rectify image sequence of camera %d first!\n', pp);
            return;
        end;
        % output
        active_imgviews{pp} = false(1,n_frame);
        flag = input('Set this camera view active or not? ([]=yes, other=no) ','s');
        if isempty(flag),
            active_imgviews{pp}(ind_active) = 1;
        end;
        imgdir_cell{pp} = imgdir;
        imgbase_cell{pp} = imgbase;
        imgfmt_cell{pp} = imgfmt;
        nframe_slot(pp) = n_frame;
        diff_imnum{pp} = diff(frame_num);
        strnum_frameNcam{pp} = strnum_frame;
    end;
    if ~all(nframe_slot==n_frame),
        fprintf(1,'\nERROR: Image numbers of every camera view are not consistent with each other!\n');
        return;
    end;
    if n_frame>1,
        for pp = 2:n_cam,
            if ~isequal(diff_imnum{pp}, diff_imnum{1}),
                fprintf(1,'\nERROR: Image number differences of camera %d are not consistent with camera 1!\n',pp);
                return;
            end;
        end;
    end;
    if any(diff(diff_imnum{1})),
        fprintf(1,'\nWARNING: It is recommended that all image frames have the same interval!\n');
    end;
    active_imgviews = cell2mat(active_imgviews);
    % set the saving directory to the first camera view folder
    imgdir = imgdir_cell{1};
end;

save_name = [imgdir '/Extraction_result.mat'];
string_save = ['save ' save_name ' fgprefix calib_mode active_imgviews n_frame n_cam handcc Omcc Tcc'];
if exist('Omcc_error','var'),
    string_save = [string_save ' Omcc_error Tcc_error'];
end;

% find the codedir of the toolbox
codeDir = mfilename('fullpath');
flag = find(codeDir=='\' | codeDir=='/',1,'last');
if ~isempty(flag),
    codeDir = codeDir(1 : flag(end)-1);
end;
% write m file load_imgdir.m into codeDir
fid = fopen([codeDir,'/load_imgdir.m'],'wt'); % text mode
fprintf(fid,'%%%% pathname for images %%%%\n\n');
fprintf(fid,'imgdir = ''%s'';\n\n', imgdir);
if calib_mode,
    string_save = [string_save ' imgbase imgfmt strnum_frame nx ny fc cc alpha_c kc'];
    if exist('fc_error','var'),
        string_save = [string_save ' fc_error cc_error alpha_c_error kc_error'];
    end;
else
    fprintf(fid,'imgdir_cell = cell(1,%d);\n', n_cam);
    for pp = 1:n_cam,
        fprintf(fid,'imgdir_cell{%d} = ''%s'';\n', pp,imgdir_cell{pp});
    end;
    string_save = [string_save ' imgbase_cell imgfmt_cell strnum_frameNcam' ...
        ' imsize fc_mat cc_mat alpha_vec kc_mat'];
    if exist('fc_mat_error','var'),
        string_save = [string_save ' fc_mat_error cc_mat_error alpha_vec_error kc_mat_error'];
    end;
end;
eval(string_save);
fclose(fid);
copyfile([codeDir,'/load_imgdir.m'], [imgdir,'/load_imgdir.m']);
fprintf(1,'\nDone with foreground extraction ...\n\n');
