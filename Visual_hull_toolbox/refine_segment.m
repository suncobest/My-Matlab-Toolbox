% track and segment members given initial segmentation
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

save_name = [imgdir '/segment_variables.mat'];
if exist(save_name,'file')==2,
    fprintf(1,'\nLoading segmentation variable file ''segment_variables.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\nSegmentation variables ''segment_variables.mat'' not found!\n');
    return;
end;

for ii=1:n_unit,
    nc = sprintf(['%0' ndigit 'd'],ii);
    save_name = [imgdir '/segment_part' nc '.mat'];
    if exist(save_name,'file')~=2,
        fprintf(1,'\nERROR: cannot find data ''%s''!\n',save_name);
        return;
    end;
end;

filepath = [imgdir '/segment'];
if exist(filepath,'dir')==7,
    save_name = [filepath '/segment_variables.mat'];
    if exist(save_name,'file')~=2,
        fprintf(1,'\nOriginal segmentation variables not found!\n');
        return;
    end;
    for ii=1:n_unit,
        nc = sprintf(['%0' ndigit 'd'],ii);
        save_name = [filepath '/segment_part' nc '.mat'];
        if exist(save_name,'file')~=2,
            fprintf(1,'\nERROR: cannot find data ''%s''!\n',save_name);
            return;
        end;
    end;
    fprintf(1,'\nOriginal segmentation data detected!\n');
    flag = input('Refine from last segmentation or original one? ([]=last, other=original) ','s');
    if ~isempty(flag),
        save_name = 'segment_variables.mat';
        copyfile([filepath '/' save_name],[imgdir '/' save_name]);
        load([imgdir '/' save_name]);
        for ii=1:n_unit,
            nc = sprintf(['%0' ndigit 'd'],ii);
            save_name = ['segment_part' nc '.mat'];
            copyfile([filepath '/' save_name],[imgdir '/' save_name]);
        end;
    end;
else
    % copy orignial data to directory 'segment'
    mkdir(imgdir,'segment');
    save_name = 'segment_variables.mat';
    copyfile([imgdir '/' save_name],[filepath '/' save_name]);
    for ii=1:n_unit,
        nc = sprintf(['%0' ndigit 'd'],ii);
        save_name = ['segment_part' nc '.mat'];
        copyfile([imgdir '/' save_name],[filepath '/' save_name]);
    end;
end;

if ~exist('show_ellipsoid','var'),
    show_ellipsoid = input('Show the principal ellipsoid for every parts or not? ([]=yes, other=no) ','s');
    show_ellipsoid = isempty(show_ellipsoid);
end;
if ~exist('nlen_vec','var'),
    nlen_vec = input('Length factor of principle axis vector to show: ([]=1.5 times semi_axes) ');
    if isempty(nlen_vec),
        nlen_vec = 1.5;
    end;
end;
MaxIter = input('Maximum iteration number for the ellipsoid_fitting function: ([]=10) ');
if isempty(MaxIter),
    MaxIter = 10;
end;
fracw = input('Fraction number for wings to estimate chord direction: ([]=0.1) ');
if isempty(fracw),
    fracw = 0.1;
end;
fracb = input('Fraction number for body to estimate heading direction: ([]=0.2) ');
if isempty(fracb),
    fracb = 0.2;
end;
flag = 1;
fprintf(1,'\nTracking function: 1 ellipsoid_tracking!\n');
while flag,
    func_flag = input('Which function to call: (1=ellipsoid_tracking, 2=ellipsoid_rotation, 3=ellipsoid_rotation2) ');
    if length(func_flag)~=1 || all(func_flag~=[1 2 3]),
        fprintf('Unexpected input! Please enter again!\n')
        continue;
    end;
    flag = 0;
end;
SW = input('Update body orientation or not? ([]=yes, other=no) ','s');
SW = isempty(SW);

Nb = 20;
npts = (Nb+1)^2;
[Xe, Ye, Ze] = sphere(Nb);
XYZs = [Xe(:)'; Ye(:)'; Ze(:)'];

ind_active = find(active_images);
N_active = length(ind_active);
assert(isequal(diff(ind_active), ones(1,N_active-1)),'All active images must be continuous to extract kinematics!');
% smoothing kinematics
Npc = 10;
nr = floor(FramePC/Npc);
fmask = fspecial('gaussian',[1,2*nr+1],nr/2);
center = reshape(imfilter(reshape(member_center,3*(nparts+1),[]),fmask),3,nparts+1,[]);
member_center(:,:,nr+1:end-nr) = center(:,:,nr+1:end-nr);
for ii = 1:nparts+1,
    Qt = quatsmooth(trans_quat_axis(reshape(member_omega(:,ii,:),3,[])),fmask);
    member_omega(:,ii,:) = reshape(trans_quat_axis(Qt),3,1,[]);
end;
% average size and main axis direction of all members
axes_mean = sum(member_axes(:,:,ind_active),3)/N_active;
axis_ux = NaN(3,nparts+1,N_active);
for kk=ind_active,
    for ii=1:nparts+1,
        R = rodrigues(member_omega(:,ii,kk));
        axis_ux(:,ii,kk) = R(:,1);
    end;
end;
meanfoot = member_center(:,:,ind_active) - repmat(axes_mean(1,:),[3,1,N_active]).*axis_ux;
body_y = reshape(meanfoot(:,2,:)-meanfoot(:,3,:), 3,[]);
body_y = body_y./repmat(sqrt(sum(body_y.*body_y)),3,1);    % axis2 of body
body_x = reshape(axis_ux(:,1,:), 3,[]);   % axis1 of body
% using middle downstroke to denoise body's lateral direction
Nb = length(ind_maxspan);
%  frame numbers of middle downstroke and upstroke
a = body_x(:,ind_maxspan);
b = body_y(:,ind_maxspan);
b = b-repmat(sum(b.*a),3,1).*a;
b = b./repmat(sqrt(sum(b.*b)),3,1);
ind = ind_maxspan;
omega = zeros(3,Nb);   % body attitudes with maximum wing span
for ii = 1:Nb,
    omega(:,ii) = rodrigues([a(:,ii), b(:,ii), cross(a(:,ii),b(:,ii))]);
end;
% quaternion interpolation (squad function) with clamped orientation
if ind(1)>ind_active(1),
    ind = [ind_active(1),ind];
    omega = [omega(:,1),omega];
end;
if ind(end)<ind_active(end),
    ind = [ind, ind_active(end)];
    omega = [omega, omega(:,end)];
end;
% update the body orientation
omega = trans_quat_axis(squad(ind, trans_quat_axis(omega),ind_active));
member_omega(:,1,ind_active) = omega;

% compute wings' center, tip, root position and orientation with respect to body coordinate system
wcenter = member_center(:,2:3,ind_active)-repmat(member_center(:,1,ind_active),[1, 2]);
XX = [wcenter,axis_ux(:,2:3,:)];
for kk = ind_active,
    XX(:,:,kk) = rodrigues(omega(:,kk))'*XX(:,:,kk);
end;
wcenter = reshape(XX(:,1:2,:),3,[]);
wing_x = reshape(XX(:,3:4,:),3,[]);
% compute wing roots
wroot(:,1) = lines_joint(wcenter(:,1:2:end), wing_x(:,1:2:end));
wroot(:,2) = lines_joint(wcenter(:,2:2:end), wing_x(:,2:2:end));
wroot(2,2) = -wroot(2,2);
wroot = sum(wroot,2)/2*[1,1];
wroot(2,2) = -wroot(2,2);

Nstart = 1;
kk = Nstart;
if ~active_images(kk),
    fprintf(1,'\nProgram need the start frame to be active!\n');
    return;
end;
count = ceil(kk/nfpu);
base = (count-1)*nfpu;
bk = kk-base;
nc = sprintf(['%0' ndigit 'd'],count);
load([imgdir '/segment_part' nc '.mat']);
XX = X_cell{bk};
seg_id = seg_cell{bk};

% compute body axes
center = member_center(:,1,kk);
vect = rodrigues(member_omega(:,1,kk));
indk = (seg_id==1) | (seg_id==0);
XXk = XX(:, indk);
% resize ellipsoid
re_edit = 1;
while re_edit,
    fprintf(1,'\nResize section ellipsoid of part %d...\n',1);
    Nb = input('Resize factor for the semi axes: ([]=1)');
    n = length(Nb);
    if n==0 || all(Nb==1),
        fprintf(1,'Do not resize at all!\n');
        semiax = axes_mean(:,1);
    elseif (n~=1 && n~=3) || any(Nb<=0),
        fprintf(1,'Unexpected input! Please enter again!\n');
        continue;
    else
        semiax = axes_mean(:,1).*Nb(:);
    end;
    [ct, ax, vec, ind] = ellipsoid_fitting(XXk,center,semiax,vect, MaxIter, SW);
    idx = ~ind;
    figure(3);
    plot3(XXk(1,ind),XXk(3,ind),-XXk(2,ind), '.', XXk(1,idx),XXk(3,idx),-XXk(2,idx), '.');
    hold on;
    % generate ellipsoid from unit sphere
    XYZe = vec*diag(semiax)*XYZs+ct(:,ones(1,npts));
    Xe(:) = XYZe(1,:);
    Ye(:) = XYZe(2,:);
    Ze(:) = XYZe(3,:);
    surf(Xe,Ze,-Ye,'FaceColor',palette(1,:), 'EdgeColor',palette(1,:), 'FaceAlpha',0.05);
    hold off;
    view(az,el);
    fprintf(1,'Check if the principal ellipsoid need to be resized:\n');
    figure(3); axis equal off;
    re_edit = input('Need to resize the ellipsoid or not? ([]=yes, other=no) ','s');
    re_edit = isempty(re_edit);
end;
temp = uint8(false(size(ind)));     % idx=~ind are treated as outliners (seg_id == 0)
temp(ind) = 1;
seg_id(indk) = temp;
member_center(:,1,kk) = ct;
member_axes(:,1,kk) = ax;
axes_mean(:,1) = semiax;    % update segment ellipsoid

% refine the heading direction
vc =  vec(:,1);
ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*(1-fracb);
XXk = XXk(:,ind);
Nb = size(XXk,2);
idx = vc'*(XXk-ctt(:,1)*ones(1,Nb))>=0;
ind = vc'*(XXk-ctt(:,2)*ones(1,Nb))<=0;
vc = mean(XXk(:,idx),2)-mean(XXk(:,ind),2);
% Schmidt orthogonalization
vc = vc/norm(vc);
vec(:,1) = vc;
body_x = vc;
vc = vec(:,2)-vec(:,2)'*vc*vc;
vec(:,2) = vc/norm(vc);
vec(:,3) = cross(vec(:,1), vec(:,2));
member_omega(:,1,kk) = rodrigues(vec);

% compute wing axes
root = vect*wroot+ct(:,[1 1]);
for ii=2:nparts+1,
    center = member_center(:,ii,kk);
    vect = rodrigues(member_omega(:,ii,kk));
    indk = (seg_id==ii) | (seg_id==0);
    XXk = XX(:, indk);
    % resize ellipsoid
    re_edit = 1;
    while re_edit,
        fprintf(1,'\nResize section ellipsoid of part %d...\n',ii);
        Nb = input('Resize factor for the semi axes: ([]=1)');
        n = length(Nb);
        if n==0 || all(Nb==1),
            fprintf(1,'Do not resize at all!\n');
            semiax = axes_mean(:,ii);
        elseif (n~=1 && n~=3) || any(Nb<=0),
            fprintf(1,'Unexpected input! Please enter again!\n');
            continue;
        else
            semiax = axes_mean(:,ii).*Nb(:);
        end;
        switch func_flag,
            case 1,
                [ct, ax, vec, ind] = ellipsoid_tracking(XXk,root(:,ii-1),semiax,vect,MaxIter, fracw);
            case 2,
                [ct, ax, vec, ind] = ellipsoid_rotation(XXk,root(:,ii-1),semiax,vect,MaxIter);
                % project wing points between two planes of a faction of span distance
                % on to the plane perpendicular to span axis and through the center
                Xk = XXk(:,ind);
                Nb = size(Xk,2);
                vc = vec(:,1);
                ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*fracw;
                idx = vc'*(Xk-ctt(:,1)*ones(1,Nb))<=0 & vc'*(Xk-ctt(:,2)*ones(1,Nb))>=0;
                Xk = Xk(:,idx);
                Xk = Xk-ct(:,ones(1,size(Xk,2)));
                Xk = Xk-vc*(vc'*Xk);
                vc = vect(:,2)-vc'*vect(:,2)*vc;
                [~, vc] = max_projection_vector(Xk,vc);
                vec(:,2) = vc/norm(vc);
                vec(:,3) = cross(vec(:,1), vec(:,2));
            case 3,
                [ct, ax, vec, ind] = ellipsoid_rotation2(XXk,root(:,ii-1),semiax,vect,MaxIter, fracw);
            otherwise,
                error('Unexpected value for function index!');
        end;
        idx = ~ind;
        figure(3);
        plot3(XXk(1,ind),XXk(3,ind),-XXk(2,ind), '.', XXk(1,idx),XXk(3,idx),-XXk(2,idx), '.');
        hold on;
        % generate ellipsoid from unit sphere
        XYZe = vec*diag(semiax)*XYZs+ct(:,ones(1,npts));
        Xe(:) = XYZe(1,:);
        Ye(:) = XYZe(2,:);
        Ze(:) = XYZe(3,:);
        surf(Xe,Ze,-Ye,'FaceColor',palette(1,:), 'EdgeColor',palette(1,:), 'FaceAlpha',0.05);
        hold off;
        view(az,el);
        fprintf(1,'Check if the principal ellipsoid need to be resized:\n');
        figure(3); axis equal off;
        re_edit = input('Need to resize the ellipsoid or not? ([]=yes, other=no) ','s');
        re_edit = isempty(re_edit);
    end;
    temp = uint8(false(size(ind)));     % idx=~ind are treated as outliners (seg_id == 0)
    temp(ind) = ii;
    seg_id(indk) = temp;
    member_center(:,ii,kk) = ct;
    member_axes(:,ii,kk) = ax;
    axes_mean(:,ii) = semiax;    % update segment ellipsoid
    member_omega(:,ii,kk) = rodrigues(vec);
end;
% save segment index
seg_cell{bk} = seg_id;
% compute lateral direction of body (from left to right wing center)
body_y = member_center(:,2,kk)-member_center(:,3,kk);
body_y = body_y-body_y'*body_x*body_x;     % orthogonalization
body_y = body_y/norm(body_y);   %  axis2 of body
% update the body orientation
member_omega(:,1,kk) = rodrigues([body_x, body_y, cross(body_x,body_y)]);
axes_mean(:,2:3) = sum(axes_mean(:,2:3),2)/2*[1,1];

plot_string = 'plot3(';
for k=1:nparts+1,
    plot_string = [plot_string 'XX(1,seg_id==' num2str(k) '),XX(3,seg_id==' ...
        num2str(k) '),-XX(2,seg_id==' num2str(k) '),''.'','];
end;
plot_string = [plot_string 'XX(1,seg_id==0),XX(3,seg_id==0),-XX(2,seg_id==0),''.'')'];
figure(3);
eval(plot_string); hold on;
for ii = 1:nparts+1,
    % plot the three semi axes
    ct = member_center(:,ii,kk);
    ax = axes_mean(:,ii);
    vc = rodrigues(member_omega(:,ii,kk));
    ctt = ct(:,ones(1,3));
    xyz = ctt+vc*diag(ax)*nlen_vec;
    for i=1:3,
        plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
    end;
    if show_ellipsoid,
        % generate ellipsoid from unit sphere
        XYZe = vc*diag(ax)*XYZs+ct(:,ones(1,npts));
        Xe(:) = XYZe(1,:);
        Ye(:) = XYZe(2,:);
        Ze(:) = XYZe(3,:);
        surf(Xe,Ze,-Ye,'FaceColor',palette(ii,:), 'EdgeColor',palette(ii,:), 'FaceAlpha',0.05);
    end;
end;
axis equal vis3d off;
hold off;
view(az,el);
set(gcf,'color',[1 1 1]*0.7);
title(['Segment result of frame: ', num2str(kk)]);
drawnow;
if kk == base+nfpu || kk == n_frame,
    fprintf(1,'Saving segmentation result ''segment_part%s.mat''.\n',nc);
    save_name = [imgdir '/segment_part' nc '.mat'];
    eval(['save ' save_name ' X_cell seg_cell']);
end;

if Nstart==n_frame,
    return;
end;

% initial speed of translation and quaternion rotaion
for kk = Nstart+1:n_frame,
    if kk == base+nfpu+1,
        base = base+nfpu;
        count = count+1;
        nc = sprintf(['%0' ndigit 'd'],count);
        % load visual hull data
        load([imgdir '/segment_part' nc '.mat']);
    end;
    if active_images(kk),
        bk = kk-base;
        XX = X_cell{bk};
        seg_id = seg_cell{bk};
        
        % recompute body axes
        center = member_center(:,1,kk);
        vect = rodrigues(member_omega(:,1,kk));
        indk = (seg_id==1) | (seg_id==0);
        XXk = XX(:, indk);
        [ct, ax, vec, ind] = ellipsoid_fitting(XXk,center,axes_mean(:,1),vect, MaxIter, SW);
        temp = uint8(false(size(ind)));     % idx=~ind are treated as outliners (seg_id == 0)
        temp(ind) = 1;
        seg_id(indk) = temp;
        member_center(:,1,kk) = ct;
        member_axes(:,1,kk) = ax;
        
        % refine the heading direction
        vc =  vec(:,1);
        ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*(1-fracb);
        XXk = XXk(:,ind);
        Nb = size(XXk,2);
        idx = vc'*(XXk-ctt(:,1)*ones(1,Nb))>=0;
        ind = vc'*(XXk-ctt(:,2)*ones(1,Nb))<=0;
        vc = mean(XXk(:,idx),2)-mean(XXk(:,ind),2);
        body_x = vc/norm(vc);
        
        % compute wing axes
        root = vect*wroot+ct(:,[1 1]);
        for ii=2:nparts+1,
            vect = rodrigues(member_omega(:,ii,kk));
            indk = (seg_id==ii) | (seg_id==0);
            XXk = XX(:, indk);
            switch func_flag,
                case 1,
                    [ct, ax, vec, ind] = ellipsoid_tracking(XXk,root(:,ii-1),axes_mean(:,ii),vect,MaxIter, fracw);
                case 2,
                    [ct, ax, vec, ind] = ellipsoid_rotation(XXk,root(:,ii-1),axes_mean(:,ii),vect,MaxIter);
                    % project wing points between two planes of a faction of span distance
                    % on to the plane perpendicular to span axis and through the center
                    XXk = XXk(:,ind);
                    Nb = size(XXk,2);
                    vc = vec(:,1);
                    ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*fracw;
                    idx = vc'*(XXk-ctt(:,1)*ones(1,Nb))<=0 & vc'*(XXk-ctt(:,2)*ones(1,Nb))>=0;
                    XXk = XXk(:,idx);
                    XXk = XXk-ct(:,ones(1,size(XXk,2)));
                    XXk = XXk-vc*(vc'*XXk);
                    vect = vect(:,2)-vc'*vect(:,2)*vc;
                    [~, vc] = max_projection_vector(XXk,vect);
                    vec(:,2) = vc/norm(vc);
                    vec(:,3) = cross(vec(:,1), vec(:,2));
                case 3,
                    [ct, ax, vec, ind] = ellipsoid_rotation2(XXk,root(:,ii-1),axes_mean(:,ii),vect,MaxIter, fracw);
                otherwise,
                    error('Unexpected value for function index!');
            end;
            temp = uint8(false(size(ind)));     % idx=~ind are treated as outliners (seg_id == 0)
            temp(ind) = ii;
            seg_id(indk) = temp;
            member_center(:,ii,kk) = ct;
            member_axes(:,ii,kk) = ax;
            member_omega(:,ii,kk) = rodrigues(vec);
        end;
        seg_cell{bk} = seg_id;
        % compute lateral direction of body (from left to right wing center)
        body_y = member_center(:,2,kk)-member_center(:,3,kk);
        body_y = body_y-body_y'*body_x*body_x;     % orthogonalization
        body_y = body_y/norm(body_y);   %  axis2 of body
        % update the body orientation
        member_omega(:,1,kk) = rodrigues([body_x, body_y, cross(body_x,body_y)]);
        
        % plot segment result
        plot_string = 'plot3(';
        for k=1:nparts+1,
            plot_string = [plot_string 'XX(1,seg_id==' num2str(k) '),XX(3,seg_id==' ...
                num2str(k) '),-XX(2,seg_id==' num2str(k) '),''.'','];
        end;
        plot_string = [plot_string 'XX(1,seg_id==0),XX(3,seg_id==0),-XX(2,seg_id==0),''.'')'];
        figure(3);
        eval(plot_string); hold on;
        for ii = 1:nparts+1,
            % plot the three semi axes
            ct = member_center(:,ii,kk);
            ax = axes_mean(:,ii);
            vc = rodrigues(member_omega(:,ii,kk));
            ctt = ct(:,ones(1,3));
            xyz = ctt+vc*diag(ax)*nlen_vec;
            for i=1:3,
                plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
            end;
            if show_ellipsoid,
                % generate ellipsoid from unit sphere
                XYZe = vc*diag(ax)*XYZs+ct(:,ones(1,npts));
                Xe(:) = XYZe(1,:);
                Ye(:) = XYZe(2,:);
                Ze(:) = XYZe(3,:);
                surf(Xe,Ze,-Ye,'FaceColor',palette(ii,:), 'EdgeColor',palette(ii,:), 'FaceAlpha',0.05);
            end;
        end;
        axis equal vis3d off;
        hold off;
        view(az,el);
        set(gcf,'color',[1 1 1]*0.7);
        title(['Segment result of frame: ', num2str(kk)]);
        drawnow;
    end;
    if kk == base+nfpu || kk == n_frame,
        fprintf(1,'Saving segmentation result ''segment_part%s.mat''.\n',nc);
        save_name = [imgdir '/segment_part' nc '.mat'];
        eval(['save ' save_name ' X_cell seg_cell']);
    end;
end;

% update segment variables
fprintf(1,'\nDone with segmentation refinement...\nUpdate segment variables in ''segment_variables.mat''.\n');
save_name = [imgdir '/segment_variables.mat'];
string_save = ['save ' save_name ' n_frame ind_active active_images ndigit palette color_cell kparts nparts '...
    'member_center member_axes member_omega nfpu n_unit FramePS cyclePS FramePC Nstroke nfpc ind_down '...
    'ind_up ind_maxdown ind_maxup ind_maxspan wroot Qbw dfpc az el'];
eval(string_save);
fprintf(1,'done...\nNow you can extract the kinematics...\n');
