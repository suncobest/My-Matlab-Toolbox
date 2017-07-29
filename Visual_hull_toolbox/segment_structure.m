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

save_name = [imgdir '/visualhull_environment.mat'];
if exist(save_name, 'file')==2,
    fprintf(1,'\nLoading environment variable file ''visualhull_environment.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\nERROR: ''visualhull_environment.mat'' not found!\n');
    return;
end;

nc = sprintf(['%0' ndigit 'd'],n_unit);
save_name = [imgdir '/segment_variables.mat'];
if exist(save_name,'file')==2,
    save_name = [imgdir '/segment_part' nc '.mat'];
    if exist(save_name, 'file')==2,
        fprintf(1,'\nOld segmentation file detected ...\n');
        flag = input('Overwrite result of segmentation or not? ([]=no, other=yes) ','s');
        if isempty(flag),
            fprintf(1,'\nSegmentation process terminated!\n');
            return;
        end;
    end;
end;

% check 3D points data
nc = sprintf(['%0' ndigit 'd'], 1);
save_name = [imgdir '/3D_points_part' nc '.mat'];
if exist(save_name, 'file')==2,
    load(save_name);
else
    fprintf(1,'\nERROR: cannot find data ''%s''!\n',save_name);
    return;
end;

show_ellipsoid = input('Show the principal ellipsoid for every parts or not? ([]=yes, other=no) ','s');
show_ellipsoid = isempty(show_ellipsoid);
nlen_vec = input('Length factor of principle axis vector to show: ([]=1.5 times semi_axes) ');
if isempty(nlen_vec),
    nlen_vec = 1.5;
end;

% mesh of unit sphere
Nb = 20;
npts = (Nb+1)^2;
[Xe, Ye, Ze] = sphere(Nb);
XYZs = [Xe(:)'; Ye(:)'; Ze(:)'];
% colors for points
palette = lines(7);     % lines returns an M-by-3 matrix containing a "ColorOrder" colormap
color_cell = {'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'dark'};

kk = 1;
% load 3D points cloud of frame 1
XX = X_cell{kk};
seg_cell = cell(size(X_cell));
if isempty(XX),
    fprintf(1,'\nERROR: No cuboids found on border of frame %d ...\n', kk);
    return;
end;

% segment test process
fprintf(1,'Initialization of points cloud segmentation:\nFirst segment the main body from points cloud...\n');
[center, semiax, vec, ind] = principal_ellipsoid(XX);
idx = ~ind;
figure(3); hold on;
plot3(XX(1,ind),XX(3,ind),-XX(2,ind),'.', 'color', palette(1,:));
plot3(XX(1,idx),XX(3,idx),-XX(2,idx),'.', 'color', palette(2,:));
% plot the three semi axes
ctt = center(:,ones(1,3));
xyz = ctt+vec*diag(semiax)*nlen_vec;
for i=1:3,
    plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
end;
axis equal vis3d off;
% view axis for 3D space
re_edit = 1;
while re_edit,
    temp = input('Set looking direction for 3D view: ([azimuth, elevation]) ');
    if length(temp)~=2,
        fprintf(1,'Unexpected input! Please enter again!\n');
        continue;
    end;
    az = temp(1);
    el = temp(2);
    figure(3);
    view(az,el);
    re_edit = input('Need to reset looking direction or not? ([]=yes, other=no) ','s');
    re_edit = isempty(re_edit);
end;

% adjust the orientation of body axes
fprintf(1, 'Identify the orientation of body axes: 1=roll axis; 2=pitch axis; 3=yaw axis;\n');
fprintf(1,'\nCheck the color and direction of roll (longitudinal) axis:\n');
figure(3);
ax_id = 1:3;
re_edit = 1;
while re_edit,
    idk = input('What is the color of axis 1? (1=blue; 2=Green; 3=red;) ');
    if ~isscalar(idk) || all(idk~=1:3),
        fprintf(1,'Unexpected input! Please enter again!\n');
        continue;
    end;
    re_edit = 0;
end;
ax_id(1) = idk;
flag = input('Does the roll axis point to head or tail? ([]=head, other=tail) ','s');
if ~isempty(flag),
    vec(:,idk) = -vec(:,idk);
end;

fprintf(1,'\nCheck the color and direction of pitch (lateral) axis:\n');
re_edit = 1;
while re_edit,
    idk = input('What is the color of axis 2? (1=blue; 2=Green; 3=red;) ');
    if ~isscalar(idk) || all(idk~=1:3) || idk==ax_id(1),
        fprintf(1,'Unexpected input! Please enter again!\n');
        continue;
    end;
    re_edit = 0;
end;
ax_id(2) = idk;
flag = input('Does the pitch axis point to right or left? ([]=right, other=left) ','s');
if ~isempty(flag),
    vec(:,idk) = -vec(:,idk);
end;
ax_id(3) = 6-sum(ax_id(1:2));
% rectify axes to right handed
vec = vec(:,ax_id);
vec(:,3) = cross(vec(:,1), vec(:,2));
semiax = semiax(ax_id);

% generate ellipsoid from unit sphere
XYZe = vec*diag(semiax)*XYZs+center(:,ones(1,npts));
Xe(:) = XYZe(1,:);
Ye(:) = XYZe(2,:);
Ze(:) = XYZe(3,:);
figure(3);
surf(Xe,Ze,-Ye,'FaceColor',palette(1,:), 'EdgeColor',palette(1,:), 'FaceAlpha',0.05);
hold off;
drawnow;

% resize the ellipsoid
fprintf(1,'\nCheck if the principal ellipsoid need to be resized:\n');
re_edit = input('Need to resize the ellipsoid or not? ([]=yes, other=no) ','s');
re_edit = isempty(re_edit);
while re_edit,
    Nb = input('Resize factor for the semi axes: ([]=1)');
    n = length(Nb);
    if n==0 || all(Nb==1),
        ax = semiax;
        fprintf(1,'Do not resize at all!\n');
    elseif (n~=1 && n~=3) || any(Nb<=0),
        fprintf(1,'Unexpected input! Please enter again!\n');
        continue;
    else
        ax = semiax.*Nb(:);
    end;
    [ct, ax, vc, ind] = principal_ellipsoid(XX,center,ax,vec);
    idx = ~ind;
    figure(3);
    plot3(XX(1,ind),XX(3,ind),-XX(2,ind),'.', 'color', palette(1,:)); hold on;
    plot3(XX(1,idx),XX(3,idx),-XX(2,idx),'.', 'color', palette(2,:));
    ctt = ct(:,ones(1,3));
    xyz = ctt+vc*diag(ax)*nlen_vec;
    for i=1:3,
        plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
    end;
    axis equal vis3d off;
    % generate ellipsoid from unit sphere
    XYZe = vc*diag(ax)*XYZs+ct(:,ones(1,npts));
    Xe(:) = XYZe(1,:);
    Ye(:) = XYZe(2,:);
    Ze(:) = XYZe(3,:);
    surf(Xe,Ze,-Ye,'FaceColor',palette(1,:), 'EdgeColor',palette(1,:), 'FaceAlpha',0.05);
    hold off;
    axis equal vis3d off;
    view(az,el);
    set(gcf,'color',[1 1 1]*0.7);
    drawnow;
    fprintf(1,'\nCheck if the principal ellipsoid need to be resized:\n');
    re_edit = input('Need to resize the ellipsoid or not? ([]=yes, other=no) ','s');
    re_edit = isempty(re_edit);
    if ~re_edit,
        center = ct;
        semiax = ax;
        vec = vc;
    end;
end;
% segment index of all points
seg_id = uint8(ind);

% segment the rest points into "kparts" parts, and choose "nparts" valid parts from them.
fprintf(1, 'Second using kmeans function to cluster rest points besides the main body...\n');
XXk = XX(:,idx);
figure(3);
plot3(XXk(1,:),XXk(3,:),-XXk(2,:), '.');
axis equal vis3d off;
view(az,el);
fprintf(1,'\nCheck rest points to be clustered:\n');
re_edit = 1;
while re_edit,
    flag = 1;
    while flag,
        kparts = input('How many parts to be segmented? ([]=3) ');
        if isempty(kparts),
            kparts = 3;
        else
            kparts = round(kparts);
            if kparts<1 || kparts>6,
                fprintf(1,'\nUnexpected input! Please enter a number in the region of 1~6!\n');
                continue;
            end;
        end;
        flag = 0;
    end;
    [indk,ct] = kmeans(XXk',kparts);
    Nb = zeros(1,kparts);
    for ii=1:kparts,
        Nb(ii) = sum(indk==ii);
    end;
    % indices of kparts in descend order of points amount
    [~,ind] = sort(Nb,'descend');
    plot_string = 'plot3(';
    for ii=1:kparts,
        k = ind(ii);
        plot_string = [plot_string 'XXk(1,indk==' num2str(k) '),XXk(3,indk==' ...
            num2str(k) '),-XXk(2,indk==' num2str(k) '),''.'','];
    end;
    plot_string = [plot_string(1:end-1) ');'];
    figure(3);
    eval(plot_string);
    view(az,el);
    axis equal vis3d off;
    fprintf(1,'\nCheck if the number of segments need to be rectified:\n');
    re_edit = input('Need to alter the number of subdivision parts or not? ([]=yes, other=no) ','s');
    re_edit = isempty(re_edit);
end;
ct = ct';
ct = ct(:,ind);  % permute parts center

% identify every valid part
fprintf(1,'\nIdentify every valid parts of the clustering results in order,\n');
flag = 1;
while flag,
    nparts = input(['Number of valid parts: ([]=2, no more than ' num2str(kparts) ') ']);
    if isempty(nparts),
        nparts = 2;
    else
        nparts = round(nparts);
        if ~isscalar(nparts) || nparts<1 || nparts>kparts,
            fprintf(1,'\nUnexpected input! Please enter again!\n');
            continue;
        end;
    end;
    flag = 0;
end;

% generate index of valid parts
flag = 1;
while flag,
    fprintf(1,'\nColor code of %d segments:\n',kparts);
    for ii=1:kparts,
        fprintf(1,[num2str(ii), '=' color_cell{ii} '; ']);
    end;
    fprintf(1,'\nPlease enter color code of %d valid parts in order: (make sure right part before left)\n',nparts);
    idk = input(['Color code vector ([]=[' num2str(1:nparts) ']): ']);
    if isempty(idk),
        idk = 1:nparts;
        flag = 0;
    else
        idk = round(idk);
        flag = length(idk)~=nparts || any(idk<1) || any(idk>kparts);
        % check idk if any color number is repeated
        if ~flag,
            for ii=1:kparts,
                flag = sum(idk==ii)>1;
                if flag, break; end;
            end;
        end;
    end;
    if flag,
        fprintf(1,'\nUnexpected input! Please enter again!\n');
    else
        fprintf(1,'\nCheck color sequence of valid parts:\n');
        for ii=1:nparts,
            fprintf(1,['Part ' num2str(ii) ': ' color_cell{idk(ii)} ';\n']);
        end;
        figure(3);
        fprintf(1,'\nCheck if the sequence of valid segments need further amendment:\n');
        flag = input('Need to alter sequence of valid parts or not? ([]=yes, other=no) ','s');
        flag = isempty(flag);
    end;
end;
% segment index for every part: 1~nparts
temp = uint8(false(1,size(XXk,2)));
for ii=1:nparts,
    temp(indk==ind(idk(ii))) = ii;
end;
% 1 are assumed as main body, nparts: 2~1+nparts
seg_id(idx) = temp+1;
% save segment index
seg_cell{kk} = seg_id;

% initialization of kinematic parameters
member_center = NaN(3,nparts+1,n_frame);
member_axes = member_center;
member_omega = member_center;
member_center(:,1,kk) = center;
member_axes(:,1,kk) = semiax;
member_omega(:,1,kk) = rodrigues(vec);

flag = 1;
while flag,
    FramePS = input('Shooting rate (FPS) of the image sequence: ');
    if isempty(FramePS) || FramePS<=0,
        fprintf(1,'\nUnexpected input!\n');
        continue;
    end;
    flag = 0 ;
end;

flag = 1;
Npc = 10;
while flag,
    FramePC = input('Average number of  frames per cycle: ');
    if isempty(FramePC) || FramePC<Npc,
        fprintf(1,'\nUnexpected input!\n');
        continue;
    end;
    flag = 0 ;
end;
cyclePS = FramePS/FramePC;      % frequency
Nstroke = (n_frame-1)/FramePC;
nfpc = round(FramePC);
dfpc = round(FramePC/5);

%% track and segment members given initial segmentation
root_mat = NaN(3,nparts,n_frame);
if Nstroke<1,
    Nstop = n_frame;
else
    Nstop = nfpc;
end;
frame_list = 1:Nstop;
Not_ok = true(1,Nstop);
while any(Not_ok),
    % initiate and track the 1st motion cycle
    track_list = frame_list(Not_ok);
    track_members;
    fprintf(1,'\nThe 1st motion cycle may need checking:\n');
    check_segmentation;
    if any(Not_ok),
        disp('Some frames need to be rectified!');
        flag = input('Interpolate these frames or not? ([]=yes, other=no)');
        if isempty(flag),
            track_list = frame_list(Not_ok);
            ind = frame_list(~Not_ok);
            ind = ind(ind>track_list(1)-3 & ind<track_list(end)+3);
            Qt = NaN(4,nparts,length(track_list));
            for ii=1:nparts,
                Qt(:,ii,:) = reshape(squad(ind, trans_quat_axis(reshape(member_omega(:,ii+1,ind),3,[])), track_list),4,1,[]);
            end;
            refine_body_wings;
            check_segmentation;
        end;
    end;
end;

if Nstroke<1,
    ind_down = [];
    ind_up = [];
    ind_maxdown = [];
    ind_maxup = [];
    ind_maxspan = [];
else
    % compute the timing of downstroke and upstroke
    vec = NaN(1,nfpc);
    for kk=1:nfpc,
        a = rodrigues(member_omega(:,2,kk));
        b = rodrigues(member_omega(:,3,kk));
        vec(kk) = sum((a(:,1)+b(:,1)).^2);
    end;
    [~,idk] = max(vec);
    ind_down = round(idk : FramePC : n_frame);
    % ind: index of downstroke; idx: index of upstroke
    if idk-FramePC/2>=1,
        ind_up = round(idk-FramePC/2 : FramePC : n_frame);
        idx = ind_up(1) : ind_down(1);
        ind = [1:ind_up(1)-1, ind_down(1)+1:nfpc];
    else
        ind_up = round(idk+FramePC/2 : FramePC : n_frame);
        ind = ind_down(1) : ind_up(1);
        idx = [1:ind_down(1)-1, ind_up(1)+1:nfpc];
    end;

    % compute the timing of wings reaching max span during downstroke
    [~,id] = min(vec(ind));
    ind_maxdown = round(ind(id) : FramePC : n_frame);
    % compute the timing of wings reaching max span during upstroke
    [~,id] = min(vec(idx));
    ind_maxup = round(idx(id) : FramePC : n_frame);
    if ind_maxdown(1)<ind_maxup(1),
        ind_maxspan = ind_maxdown;
    else
        ind_maxspan = ind_maxup;
    end;
    Nmspan = length(ind_maxspan);
end;

% Qbw: quaternion rotation from body to wings
Qbw = NaN(4,nparts,n_frame);
for kk = 1:Nstop,
    Qt = trans_quat_axis(member_omega(:,:,kk));
    Qt(1:3,1) = -Qt(1:3,1);
    Qbw(:,:,kk) = quatmul(Qt(:,1), Qt(:,2:end));
end;

nr = floor(FramePC/Npc);
fmask = fspecial('gaussian',[1,2*nr+1],nr/2);
if Nstop<n_frame,
    % apply periodic constraint to estimate only wing's incidence
    period_sw = 0;
    % set tracking function as ellipsoid_tracking: 1=ellipsoid_tracking, 2=ellipsoid_rotation, 3=ellipsoid_rotation2
    func_flag = 1;
    Nstart = Nstop+1;
    idk = ind_maxspan(1)+nfpc;
    Nstop = min(idk+dfpc, n_frame);
    track_list = Nstart:Nstop;
    idx = [nfpc-nr+1:nfpc,1:nfpc,1:nr+1];
    t0 = (0:nfpc)/nfpc;
    tn = mod(0:Nstop-Nstart,nfpc)/nfpc;
    for ii=1:nparts,
        Qt = quatsmooth(reshape(Qbw(:,ii,idx),4,[]),fmask,0);
        Qbw(:,ii,track_list) = reshape(squad(t0, Qt, tn),4,1,[]);
    end;
    % estimate direction by the last period
    track_body_wings;
    % update Qbw
    for ii =  track_list,
        Qt = trans_quat_axis(member_omega(:,:,ii));
        Qt(1:3,1) = -Qt(1:3,1);
        Qbw(:,:,ii) = quatmul(Qt(:,1), Qt(:,2:end));
    end;
    % update image number of maximum span
    if Nmspan>1,
        ind = idk-dfpc : Nstop;
        n = length(ind);
        vec = NaN(1,n);
        for ii=1:n,
            a = rodrigues(member_omega(:,2,ind(ii)));
            b = rodrigues(member_omega(:,3,ind(ii)));
            vec(ii) = sum((a(:,1)+b(:,1)).^2);
        end;
        [~,id] = min(vec);
        ind_maxspan(2) = ind(id);
        nk = 2;
    end;
end;

while Nstop<n_frame,
    if ind_maxspan(nk-1)-nr>=1 && ind_maxspan(nk)+nr<=Nstop,
        idx = ind_maxspan(nk-1)-nr:ind_maxspan(nk)+nr;
    else
        idx = [ind_maxspan(nk)-nr:ind_maxspan(nk)-1,ind_maxspan(nk-1):ind_maxspan(nk),...
             ind_maxspan(nk-1)+1:ind_maxspan(nk-1)+nr];
    end;
    nstep = ind_maxspan(nk)-ind_maxspan(nk-1);
    Nstart = Nstop+1;
    idk = ind_maxspan(nk)+nstep;
    Nstop = min(idk+dfpc, n_frame);
    track_list = Nstart:Nstop;
    t0 = (0:nstep)/nstep;
    tn = mod(track_list-ind_maxspan(nk),nstep)/nstep;
    for ii=1:nparts,
        Qt = quatsmooth(reshape(Qbw(:,ii,idx),4,[]),fmask,0);
        Qbw(:,ii,track_list) = reshape(squad(t0, Qt, tn),4,1,[]);
    end;
    % estimate direction by the last period
    track_body_wings;
    % update Qbw
    for ii =  track_list,
        Qt = trans_quat_axis(member_omega(:,:,ii));
        Qt(1:3,1) = -Qt(1:3,1);
        Qbw(:,:,ii) = quatmul(Qt(:,1), Qt(:,2:end));
    end;
    % update image number of maximum span
    if Nmspan>nk,
        ind = idk-dfpc : Nstop;
        n = length(ind);
        vec = NaN(1,n);
        for ii=1:n,
            a = rodrigues(member_omega(:,2,ind(ii)));
            b = rodrigues(member_omega(:,3,ind(ii)));
            vec(ii) = sum((a(:,1)+b(:,1)).^2);
        end;
        [~,id] = min(vec);
        ind_maxspan(nk+1) = ind(id);
        nk = nk+1;
    end;
end;

if ind_maxdown(1)==ind_maxspan(1),
    ind_maxdown = ind_maxspan;
    nstep = nfpc;
    Nb = length(ind_maxup);
    for ii=2 : Nb,
        id = ind_maxup(ii-1)+nstep;
        ind = id-dfpc : min(id+dfpc, n_frame);
        n = length(ind);
        vec = NaN(1,n);
        for kk=1:n,
            a = rodrigues(member_omega(:,2,ind(kk)));
            b = rodrigues(member_omega(:,3,ind(kk)));
            vec(kk) = sum((a(:,1)+b(:,1)).^2);
        end;
        [~,id] = min(vec);
        ind_maxup(ii) = ind(id);
        nstep = ind_maxup(ii)-ind_maxup(ii-1);
    end;
    ind_maxspan = reshape([ind_maxdown(1:Nb); ind_maxup],1,[]);
    if ind_maxdown(end)>ind_maxspan(end),
        ind_maxspan = [ind_maxspan ind_maxdown(end)];
    end;
else
    assert(ind_maxup(1)==ind_maxspan(1), 'Unexpected result!');
    ind_maxup = ind_maxspan;
    nstep = nfpc;
    Nb = length(ind_maxdown);
    for ii=2 : Nb,
        id = ind_maxdown(ii-1)+nstep;
        ind = id-dfpc : min(id+dfpc, n_frame);
        n = length(ind);
        vec = NaN(1,n);
        for kk=1:n,
            a = rodrigues(member_omega(:,2,ind(kk)));
            b = rodrigues(member_omega(:,3,ind(kk)));
            vec(kk) = sum((a(:,1)+b(:,1)).^2);
        end;
        [~,id] = min(vec);
        ind_maxdown(ii) = ind(id);
        nstep = ind_maxdown(ii)-ind_maxdown(ii-1);
    end;
    ind_maxspan = reshape([ind_maxup(1:Nb); ind_maxdown],1,[]);
    if ind_maxup(end)>ind_maxspan(end),
        ind_maxspan = [ind_maxspan ind_maxup(end)];
    end;
end;

% saving segment variables
fprintf(1,'\nDone with tracking process...\nSaving segment variables in ''segment_variables.mat''.\n');
save_name = [imgdir '/segment_variables.mat'];
string_save = ['save ' save_name ' n_frame ind_active active_images ndigit palette color_cell kparts nparts '...
    'member_center member_axes member_omega nfpu n_unit FramePS cyclePS FramePC Nstroke nfpc ind_down '...
    'ind_up ind_maxdown ind_maxup ind_maxspan wroot Qbw dfpc az el'];
eval(string_save);
fprintf(1,'done...\nNow you can refine sementation results or extract the kinematics...\n');
