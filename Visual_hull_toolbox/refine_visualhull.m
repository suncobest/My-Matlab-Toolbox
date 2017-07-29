% further refine the visual hull after initialization (See visualhull_carving)

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
    fprintf(1,'\nImage directory not found!\nYou must redirect image path in ''load_imgdir.m'' ...\n');
    return;
end;

save_name = [imgdir '/visualhull_environment.mat'];
if exist(save_name,'file')==2,
    fprintf(1,'\nLoading common data from file ''visualhull_environment.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\n''visualhull_environment.mat'' not found! You must run script ''init_visualhull'' first!\n');
    return;
end;

save_name = [imgdir '/visual_hull_level' num2str(level) '.mat'];
if exist(save_name,'file')==2,
    fprintf(1,'\nLoading file ''visual_hull_level%d.mat'' ...\n',level);
    load(save_name);
else
    fprintf(1,'\n''visual_hull_level%d.mat'' not found!\n', level);
    return;
end;

% bricksize of last level
bricksize = (boundsize./subnxyz)/(n_subdiv^(level-1));

if ~exist('saveIM','var'),
    saveIM = input('Draw and save 2D, 3D images or not? ([]=no, other=yes) ','s');
    saveIM = ~isempty(saveIM);
end;

% initial 2*2 matrix for perspective image corner (surface function)
xc = [0 1;0 1];
yc = [0 0;1 1];
zc = [1 1;1 1];
% view axis for 3D space
az = 50;
el = 45;

% format indicator of visual hull data
parted_data = 0;
% maximum number of total border points to store as one data file visual_hull_level%d.mat
if ~exist('Max_npts','var'),
    Max_npts = 1e7;   %5e7
end;

if ~exist('nfpu','var'),
    fprintf(1,'\nSome (mostly final) level will be divided into small units in case of memory deficiency!\n');
    nfpu = input('Set number of frames to save separately in one unit: ([]=50) ');
    if isempty(nfpu),
        nfpu = 50;
    end;
    n_unit = ceil(n_frame/nfpu);
end;

if calib_mode,
    if ~exist('resolution','var'),
        resolution = round(ny/height);   % dpi
    end;
    if ~exist('tpx','var') || ~exist('tpy','var'),
        tpy = ny/20;
        tpx = nx-tpy; % nx/2;
    end;
end;

for nl = level+1 : 1+n_sublevel,
    % number of bricks of last level
    npts_inhull = cellfun(@(x) size(x,2), IN_hull);
    npts_onhull = cellfun(@(x) size(x,2), ON_hull);
    % number of bricks and verts of current level
    nbricks_inhull = npts_inhull*n_subdiv^3;
    nbricks_onhull = npts_onhull*n_subdiv^3;
    %     nverts_inhull = nbricks_inhull*8;
    nverts_onhull = nbricks_onhull*8;
    sum_nverts = sum(nverts_onhull);
    
    % voxel center of current level and their cube vertices
    INmat = cell2mat(IN_hull);
    ONmat = cell2mat(ON_hull);
    if isempty(ONmat),
        fprintf(1,'\nERROR: No cuboids found on border!\n');
        return;
    end;
    % check if sum_nverts too big to save in one piece
    if sum_nverts>Max_npts,
        parted_data = 1;
        break;
    end;
    
    % Brick size of current level
    subbricksize = bricksize/n_subdiv;
    center_inhull = [];
    if ~isempty(INmat),
        center_inhull = subvoxcenter(INmat, n_subdiv, bricksize);
        %         vertices_inhull = gen_cuboids(center_inhull, subbricksize);
    end;
    center_onhull = subvoxcenter(ONmat, n_subdiv, bricksize);
    vertices_onhull = gen_cuboids(center_onhull, subbricksize);
    % process inhull points
    IN_hull = cell(1, n_frame);
    ON_hull = IN_hull;
    if ~isempty(center_inhull),
        jj2 = 0;
        for kk = ind_active,
            jj1 = jj2+1;
            jj2 = jj2 + nbricks_inhull(kk);
            IN_hull{kk} = center_inhull(:, jj1:jj2);
        end;
    end;
    
    if calib_mode,
        xx = NaN(2*n_cam, sum_nverts);
        for pp = 1:n_cam,
            om = Omcc(:,pp);
            T = Tcc(:,pp);
            hand = handcc(pp);
            x_kk = project_points_mirror2(vertices_onhull,om,T,hand,fc,cc,zeros(5,1),alpha_c);
            ii = (pp-1)*2;
            xx(ii+1:ii+2, :) = x_kk;
        end;
        
        % plot projection of 3D cuboids in every frame
        jj2 = 0; mm2 = 0;
        for kk = ind_active,
            jj1 = jj2+1;
            jj2 = jj2 + nverts_onhull(kk);
            % if a view is inactive, then all points cannot be eliminated in that view
            % default no cube in hull, all cubes on the border
            inner = true(n_cam, nbricks_onhull(kk));
            border = inner;
            for pp = 1:n_cam,
                if active_imgviews(pp,kk),
                    kth = (kk-1)*n_cam+pp;
                    ii = (pp-1)*2;
                    x_kk = xx(ii+1:ii+2, jj1:jj2)-repmat(bounding_mat(1:2,kth),1,nverts_onhull(kk))+1.5;
                    [in, on] = convhull_in_region(foreground_cell{kth}, x_kk, 8);
                    inner(pp,:) = in;
                    border(pp,:) = on | in;
                end;
            end;
            % visual hull data
            inner = all(inner,1);
            border = all(border,1) & (~inner);
            mm1 = mm2+1;
            mm2 = mm2 + nbricks_onhull(kk);
            X_kk = center_onhull(:, mm1:mm2);
            IN_hull{kk} = [IN_hull{kk}, X_kk(:, inner)];
            ON_hull{kk} = X_kk(:, border);
            
            if saveIM,
                % face index of on border cuboids
                nb = sum(border);
                subface = reshape(repmat(cface(:),1,nb)+repmat(8*(0:nb-1),24,1), 4, []);
                % draw 2D projection of border cuboids in all views
                frame_kk = false(ny,nx);
                for pp = 1:n_cam,
                    if active_imgviews(pp,kk),
                        kth = (kk-1)*n_cam+pp;
                        temp = bounding_mat(:,kth);
                        temp(1:2) = ceil(temp([2 1]));
                        frame_kk(temp(1) : temp(1)+temp(4)-1, temp(2) : temp(2)+temp(3)-1) = foreground_cell{kth};
                    end;
                end;
                figure(2); hold off;
                image(frame_kk);
                colormap(gray(2));      % colormap for binary image
                hold on;
                framenb = sprintf(['%0' ndigit 'd'],kk);
                text(tpx,tpy,['\it\fontname{Arial}\fontsize{9}\color{white}Level \color{yellow}' ...
                    num2str(nl) '\color{white} - \color{yellow}',framenb],'HorizontalAlignment','right');
                xxj = xx(:, jj1:jj2)+1;
                xxj = xxj(:, reshape(border(ones(8,1),:),1,[]));          % logical indices
                for pp = 1:n_cam,
                    if active_imgviews(pp,kk),
                        ii = (pp-1)*2;
                        x_kk = xxj(ii+1:ii+2, :);
                        patch('Faces', subface', 'Vertices', x_kk', 'FaceColor', 'b', 'EdgeColor','r','FaceAlpha',0.05);
                    end;
                end;
                % set axis position in window, position = [left, bottom, width, height]
                set(gca,'position',[0 0 1 1]);
                set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                drawnow;
                save_name = [imgdir '/level' num2str(nl) '_' framenb '.jpg'];
                print(2,'-djpeg',['-r' num2str(resolution)],save_name);
                
                % plot visual hull and perspective image in 3D space
                figure(3); hold off;
                X_kk = [IN_hull{kk},ON_hull{kk}];
                plot3(X_kk(1,:),X_kk(3,:),-X_kk(2,:),'g.');
                hc = center3d_mat(:,kk);
                frame_kk = imread(save_name);
                for pp = 1:n_cam,
                    if active_imgviews(pp,kk),
                        om = Omcc(:,pp);
                        T = Tcc(:,pp);
                        hand = handcc(pp);
                        % direction of camera frame axes in world frame: Rkk
                        % camera center in world frame: -Rkk*T
                        Rkk = rodrigues(-om);       % transpose
                        if hand~=1,
                            Rkk(3,:) = -Rkk(3,:);
                        end;
                        % positon of pespective image center
                        vect = hc+Rkk*T;
                        X_kk = hc+img_position*vect/norm(vect);
                        % position of pespective image corners (o, y, xy, x)
                        X_kk = X_kk(:,ones(1,4))+(Rkk(:,1)*[-1,-1,1,1]+Rkk(:,2)*[-1,1,1,-1])*maxbound*0.6;
                        % project the four corner on image
                        x_kk = project_points_mirror2(X_kk,om,T,hand,fc,cc,zeros(5,1),alpha_c);
                        Ikk = homography_image(frame_kk, x_kk+1, [], 0);   % image background is black
                        % transform to matlab 3D axes for surface function (o, y, x, xy)
                        X_kk = X_kk(:, [1,2,4,3]);
                        xc(:) = X_kk(1,:);
                        yc(:) = X_kk(2,:);
                        zc(:) = X_kk(3,:);
                        surface(xc, zc, -yc, Ikk,'FaceColor','texturemap','EdgeColor','none');
                    end;
                end;
                % set 3D plot
                figure(3);
                X_kk = vertices_onhull(:, jj1:jj2);
                X_kk = X_kk(:, reshape(border(ones(8,1),:),1,[]));    % border boxes
                patch('Faces', subface', 'Vertices', [X_kk(1,:); X_kk(3,:); -X_kk(2,:)]', 'FaceColor', 'b','EdgeColor','r','FaceAlpha',0.1);
                axis equal tight vis3d off;
                view(az,el);
                %  lightangle(az-90,el);  set(3, 'renderer','zbuffer');
                set(3,'color',[0 0 0], 'PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                drawnow;
                save_name = [imgdir '/3D_level' num2str(nl) '_' framenb '.jpg'];
                print(3,'-djpeg',['-r' num2str(resolution)],save_name);
            end;
        end;
        
    else
        
        xx = NaN(2*n_cam, sum_nverts);
        for pp = 1:n_cam,
            fc = fc_mat(:,pp);
            cc = cc_mat(:,pp);
            alpha_c = alpha_vec(pp);
            om = Omcc(:,pp);
            T = Tcc(:,pp);
            hand = handcc(pp);
            x_kk = project_points_mirror2(vertices_onhull,om,T,hand,fc,cc,zeros(5,1),alpha_c);
            ii = (pp-1)*2;
            xx(ii+1:ii+2, :) = x_kk;
        end;
        
        % plot projection of 3D cuboids in every frame
        jj2 = 0; mm2=0;
        for kk = ind_active,
            jj1 = jj2+1;
            jj2 = jj2 + nverts_onhull(kk);
            inner = true(n_cam, nbricks_onhull(kk));
            border = inner;
            for pp = 1:n_cam,
                if active_imgviews(pp,kk),
                    kth = (kk-1)*n_cam+pp;
                    ii = (pp-1)*2;
                    x_kk = xx(ii+1:ii+2, jj1:jj2)-repmat(bounding_mat(1:2,kth),1,nverts_onhull(kk))+1.5;
                    [in, on] = convhull_in_region(foreground_cell{kth}, x_kk, 8);
                    inner(pp,:) = in;
                    border(pp,:) = on | in;
                end;
            end;
            inner = all(inner,1);
            border = all(border,1) & (~inner);
            mm1 = mm2+1;
            mm2 = mm2 + nbricks_onhull(kk);
            X_kk = center_onhull(:, mm1:mm2);
            IN_hull{kk} = [IN_hull{kk},X_kk(:, inner)];
            ON_hull{kk} = X_kk(:, border);
            
            if saveIM,
                framenb = sprintf(['%0' ndigit 'd'],kk);
                nb = sum(border);
                subface = reshape(repmat(cface(:),1,nb)+repmat(8*(0:nb-1),24,1), 4, []);
                xxj = xx(:, jj1:jj2)+1;
                xxj = xxj(:, reshape(border(ones(8,1),:),1,[]));          % logical indices
                % plot visual hull and perspective image in 3D space
                figure(3); hold off;
                X_kk = [IN_hull{kk},ON_hull{kk}];
                plot3(X_kk(1,:),X_kk(3,:),-X_kk(2,:),'g.');
                hc = center3d_mat(:,kk);
                for pp = 1:n_cam,
                    if active_imgviews(pp,kk),
                        kth = (kk-1)*n_cam+pp;
                        nx = imsize(1,pp);
                        ny = imsize(2,pp);
                        frame_kk = false(ny,nx);
                        temp = bounding_mat(:,kth);
                        temp(1:2) = ceil(temp([2 1]));
                        frame_kk(temp(1) : temp(1)+temp(4)-1, temp(2) : temp(2)+temp(3)-1) = foreground_cell{kth};
                        % draw 2D projection of border cuboids in all camera views
                        figure(2); hold off;
                        image(frame_kk);
                        colormap(gray(2));      % colormap for binary image
                        hold on;
                        % text position. unit: points
                        tpy = ny/20;
                        tpx = nx-tpy;  % nx/2;
                        resolution = round(ny/height);   % dpi
                        text(tpx,tpy,['\it\fontname{Arial}\fontsize{9}\color{white}Level \color{yellow}' ...
                            num2str(nl) '\color{white} - Cam \color{yellow}' num2str(pp) ...
                            '\color{white} - \color{yellow}', framenb],'HorizontalAlignment','right');
                        ii = (pp-1)*2;
                        x_kk = xxj(ii+1:ii+2, :);
                        patch('Faces', subface', 'Vertices', x_kk', 'FaceColor', 'b', 'EdgeColor','r','FaceAlpha',0.05);
                        % set axis position in window, position = [left, bottom, width, height]
                        set(gca,'position',[0 0 1 1]);
                        set(2,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                        drawnow;
                        save_name = [imgdir_cell{pp} '/level' num2str(nl) '_cam' num2str(pp) '_' framenb '.jpg'];
                        print(2,'-djpeg',['-r' num2str(resolution)],save_name);
                        
                        % draw 3D
                        frame_kk = imread(save_name);
                        fc = fc_mat(:,pp);
                        cc = cc_mat(:,pp);
                        alpha_c = alpha_vec(pp);
                        om = Omcc(:,pp);
                        T = Tcc(:,pp);
                        hand = handcc(pp);
                        % direction of camera frame axes in world frame: Rkk
                        % camera center in world frame: -Rkk*T
                        Rkk = rodrigues(-om);       % transpose
                        if hand~=1,
                            Rkk(3,:) = -Rkk(3,:);
                        end;
                        % positon of pespective image center
                        vect = hc+Rkk*T;
                        X_kk = hc+img_position*vect/norm(vect);
                        % position of pespective image corners (o, y, xy, x)
                        X_kk = X_kk(:,ones(1,4))+(Rkk(:,1)*[-1,-1,1,1]+Rkk(:,2)*[-1,1,1,-1])*maxbound*0.6;
                        % project the four corner on image
                        x_kk = project_points_mirror2(X_kk,om,T,hand,fc,cc,zeros(5,1),alpha_c);
                        Ikk = homography_image(frame_kk, x_kk+1, [], 0);   % image background is black
                        % transform to matlab 3D axes for surface function (o, y, x, xy)
                        X_kk = X_kk(:, [1,2,4,3]);
                        xc(:) = X_kk(1,:);
                        yc(:) = X_kk(2,:);
                        zc(:) = X_kk(3,:);
                        figure(3);
                        surface(xc, zc, -yc, Ikk,'FaceColor','texturemap','EdgeColor','none');
                    end;
                end;
                % set 3D plot
                figure(3);
                X_kk = vertices_onhull(:, jj1:jj2);
                X_kk = X_kk(:, reshape(border(ones(8,1),:),1,[]));
                patch('Faces', subface', 'Vertices', [X_kk(1,:); X_kk(3,:); -X_kk(2,:)]', 'FaceColor', 'b','EdgeColor','r','FaceAlpha',0.1);
                axis equal tight vis3d off;
                view(az,el);
                nx = imsize(1,1);
                ny = imsize(2,1);
                resolution = round(ny/height);
                %  lightangle(az-90,el);  set(3, 'renderer','zbuffer');
                set(3,'color',[0 0 0], 'PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
                drawnow;
                % set the saving directory to the first camera view folder
                imgdir = imgdir_cell{1};
                save_name = [imgdir '/3D_level' num2str(nl) '_' framenb '.jpg'];
                print(3,'-djpeg',['-r' num2str(resolution)],save_name);
            end;
        end;
        imgdir = imgdir_cell{1};
    end;
    % update
    bricksize = subbricksize;
    level = level+1;
    save_name = [imgdir '/visual_hull_level' num2str(nl) '.mat'];
    string_save = ['save ' save_name ' IN_hull ON_hull'];
    fprintf(1,'Saving visual hull data file ''visual_hull_level%d.mat''.\n', nl);
    eval(string_save);
end;

% save separate data for further process
INmat = IN_hull;
ONmat = ON_hull;
for count = 1:n_unit,
    base = (count-1)*nfpu;
    ind = base+1 : min(base+nfpu, n_frame);
    nc = sprintf(['%0' ndigit 'd'],count);
    IN_hull = INmat(ind);
    ON_hull = ONmat(ind);
    fprintf(1,'\nSaving parted visualhull data ''level%d_part%s.mat''.\n', level,nc);
    save_name = [imgdir '/level' num2str(level) '_part' nc '.mat'];
    eval(['save ' save_name ' IN_hull ON_hull']);
end;

if parted_data,
    % continue to run refine process in memory efficient version
    refine_visuahull_part;
end;

% saving environment variables
fprintf(1,'\nUpdating environment variables in ''visualhull_environment.mat''.\n');
save_name = [imgdir '/visualhull_environment.mat'];
string_save = ['save ' save_name ' calib_mode fgprefix n_frame n_cam Tcc Omcc handcc ind_active active_images' ...
    ' active_imgviews character_len nxyz boundsize img_position maxbound bounding_mat area_mat foreground_cell' ...
    ' center2d_mat center3d_mat ex err_std ndigit cface height subnxyz level bricksize n_sublevel n_subdiv n_subdiv2 nfpu n_unit'];
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
fprintf(1,'done...\nNow refine the visual hull with final step to get 3D points cloud ...\n');

%   final refinement
compute_points_cloud;
fprintf(1,'Done!\n');
