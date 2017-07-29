%  find the codedir of the toolbox
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

save_name = [imgdir '/visualhull_environment.mat'];
if exist(save_name, 'file')==2,
    fprintf(1,'\nLoading environment variable file ''visualhull_environment.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\nERROR: ''visualhull_environment.mat'' not found!\n');
    return;
end;

fprintf(1,['\nThis scipt will show and save (2D and 3D) images of reconstruction results!\n' ...
    'Like cuboids of a exist level (0~%d) or points cloud result ...\n'],1+n_sublevel);
flag = input('Which structure do you want to check? ([]=visualhull level, other=points cloud) ','s');
flag = isempty(flag);

if flag,
    Nb = input(['Which level of cuboids structure to show? (0~' num2str(1+n_sublevel) ', []=last level) ']);
    if isempty(Nb),
        Nb = 1+n_sublevel;
    elseif Nb<0 || Nb >1+n_sublevel,
        fprintf(1,'\nERROR: input shoud be in the range of ''0~%d''!\n',1+n_sublevel);
        return;
    end;
    if Nb>0,
        temp = false(1,n_unit);
        for count = 1:n_unit,
            nc = sprintf(['%0' ndigit 'd'],count);
            save_name = [imgdir '/level' num2str(Nb) '_part' nc '.mat'];
            temp(count) = exist(save_name, 'file')==2;
        end;
        if any(temp),
            parted_data = 1;
            if ~all(temp),
                fprintf(1,'\nERROR: parted visualhull data of level %d not intact!\n',Nb);
                return;
            end;
        else
            parted_data = 0;
            save_name = [imgdir '/visual_hull_level' num2str(Nb) '.mat'];
            if exist(save_name, 'file')==2,
                fprintf(1,'\nLoading file ''visual_hull_level%d.mat'' ...\n',Nb);
                load(save_name);
            else
                fprintf(1,'\nERROR: can not find visualhull data of level %d!\n',Nb);
                return;
            end;
        end;
        % brick size of a level
        bricksize = (boundsize./subnxyz)/(n_subdiv^(Nb-1));
        fprintf(1,'\nSet switch to draw points on hull or points in hull ...\n');
        SW = input('Which points to show? ([]=points on hull, other=points in hull) ','s');
        SW = isempty(SW);
    end;
    % showing cuboids or not
    draw_faces = input('Draw cuboid faces or not? ([]=yes, other=no) ','s');
    draw_faces = isempty(draw_faces);
    if draw_faces,
        Rd = input('Render 3D faces with light or not?([]=yes, other=no) ','s');
        Rd = isempty(Rd);
    end;
else
    % check 3D points data
    for count = 1:n_unit,
        nc = sprintf(['%0' ndigit 'd'],count);
        save_name = [imgdir '/3D_points_part' nc '.mat'];
        if exist(save_name, 'file')~=2,
            fprintf(1,'\nERROR: can not find points cloud data ''%s''!\n', save_name);
            return;
        end;
    end;
end;

pim = input('Place 2D images in 3D scene or not?([]=yes, other=no) ','s');
pim = isempty(pim);

% initial 2*2 matrix for perspective image corner (surface function)
xc = [0 1;0 1];
yc = [0 0;1 1];
zc = [1 1;1 1];
% view axis for 3D space
az = 50;
el = 45;
mk = 30;

if calib_mode,
    if ~exist('resolution','var'),
        resolution = round(ny/height);   % dpi
    end;
    if ~exist('tpx','var') || ~exist('tpy','var'),
        tpy = ny/20;
        tpx = nx-tpy; % nx/2;
    end;
end;

if flag,
    if Nb==0,
        show_level0;
    else
        if parted_data,
            show_level_part;
        else
            show_visualhull_level;
        end;
    end;
else
    show_points_cloud;
end;

fprintf(1,'\nNote:\nPrefix for 2D figures: ''map_'';\nPrefix for 3D figures: ''hull_'';\n');
