% track and segment members given initial segmentation
Nstart = 1;
kk = Nstart;
count = ceil(kk/nfpu);
base = (count-1)*nfpu;
bk = kk-base;
nc = sprintf(['%0' ndigit 'd'],count);
save_name = [imgdir '/segment_part' nc '.mat'];
if exist(save_name,'file')==2,
    load(save_name);
else
    load([imgdir '/3D_points_part' nc '.mat']);
end;
XX = X_cell{bk};
seg_id = seg_cell{bk};

MaxIter = input('Maximum iteration number for the principal_ellipsoid function: ([]=50) ');
if isempty(MaxIter),
    MaxIter = 50;
end;

fracw = input('Fraction number for wings to estimate chord direction: ([]=0.1) ');
if isempty(fracw),
    fracw = 0.1;
end;
fracb = input('Fraction number for body to estimate heading direction: ([]=0.3) ');
if isempty(fracb),
    fracb = 0.3;
end;
pca_sw = input('Method to estimate wings'' chord direction: ([]=PCA, other=max_distance) ','s');
pca_sw = isempty(pca_sw);

thresh_dist = input('Distance threshold among members beside body: ([]=0.1) ');
if isempty(thresh_dist),
    thresh_dist = 0.1;
end;

% compute body axes
center = member_center(:,1,kk);
semiax = member_axes(:,1,kk);
vec = rodrigues(member_omega(:,1,kk));
indk = (seg_id==1) | (seg_id==0);
XXk = XX(:, indk);
% resize ellipsoid
re_edit = 1;
while re_edit,
    fprintf(1,'\nResize section ellipsoid of part %d...\n',1);
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
    [ct, ax, vc, ind, Nax] = principal_ellipsoid(XXk,center,ax,vec);
    idx = ~ind;
    figure(3);
    plot3(XXk(1,ind),XXk(3,ind),-XXk(2,ind), '.', XXk(1,idx),XXk(3,idx),-XXk(2,idx), '.');
    hold on;
    % generate ellipsoid from unit sphere
    XYZe = vc*diag(ax)*XYZs+ct(:,ones(1,npts));
    Xe(:) = XYZe(1,:);
    Ye(:) = XYZe(2,:);
    Ze(:) = XYZe(3,:);
    surf(Xe,Ze,-Ye,'FaceColor',palette(1,:), 'EdgeColor',palette(1,:), 'FaceAlpha',0.05);
    hold off;
    view(az,el);
    fprintf(1,'Check if the principal ellipsoid need to be resized:\nPress any key to continue...\n');
    figure(3); axis equal off;
    pause;
    re_edit = input('Need to resize the ellipsoid or not? ([]=yes, other=no) ','s');
    re_edit = isempty(re_edit);
end;
vec = vc;
temp = uint8(false(size(ind)));     % idx=~ind are treated as outliners (seg_id == 0)
temp(ind) = 1;
seg_id(indk) = temp;
member_nax(1,kk) = Nax;
member_center(:,1,kk) = ct;
member_axes(:,1,kk) = ax;

% refine the heading direction
vc =  vec(:,1);
ctt = center(:,[1,1])+vc*[1,-1]*semiax(1)*(1-fracb)/Nax;
XXk = XXk(:,ind);
Nb = size(XXk,2);
idx = vc'*(XXk-ctt(:,1)*ones(1,Nb))>=0;
ind = vc'*(XXk-ctt(:,2)*ones(1,Nb))<=0;
vc = mean(XXk(:,idx),2)-mean(XXk(:,ind),2);
% Schmidt orthogonalization
vc = vc/norm(vc);
vec(:,1) = vc;
vc = vec(:,2)-vec(:,2)'*vc*vc;
vec(:,2) = vc/norm(vc);
vec(:,3) = cross(vec(:,1), vec(:,2));
member_omega(:,1,kk) = rodrigues(vec);

indk = ~(ind | idx);
figure(3);
plot3(XXk(1,indk),XXk(3,indk),-XXk(2,indk), '.', XXk(1,ind),XXk(3,ind),-XXk(2,ind), '.', XXk(1,idx),XXk(3,idx),-XXk(2,idx), '.');
hold on;
ctt = center(:,ones(1,3));
xyz = ctt+vec*diag(semiax)*nlen_vec;
for i=1:3,
    plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
end;
view(az,el); hold off;
axis equal vis3d off;
fprintf(1,'\nCheck the definition of heading direction (blue axis):\nPress any key to continue...\n');
pause;

% compute wing axes
for ii=2:nparts+1,
    plot_string = 'plot3(';
    for k=1:nparts+1,
        plot_string = [plot_string 'XX(1,seg_id==' num2str(k) '),XX(3,seg_id==' ...
            num2str(k) '),-XX(2,seg_id==' num2str(k) '),''.'','];
    end;
    plot_string = [plot_string 'XX(1,seg_id==0),XX(3,seg_id==0),-XX(2,seg_id==0),''.'')'];
    figure(3); hold off;
    eval(plot_string);
    axis equal vis3d off;
    fprintf(1,'\nCheck shape of part %d (%s):\nPress any key to generate principal ellipsoid...\n',ii,color_cell{ii});
    pause;
    % generate ellipsoid from unit sphere
    indk = (seg_id==ii);
    XXk = XX(:, indk);
    [center, semiax, vec] = equivalent_ellipsoid(XXk);
    % re-segment the part
    indk = indk | (seg_id==0);
    XXk = XX(:, indk);
    % Generate principal ellipsoid
    re_edit = 1;
    while re_edit,
        Nb = input('Resize factor for the semi axes: ([]=1)');
        if isempty(Nb) || Nb==1,
            ax = semiax;
            fprintf(1,'Do not resize at all!\n');
        elseif Nb<=0,
            fprintf(1,'Unexpected input! Please enter again!\n');
            continue;
        else
            ax = semiax*Nb;
        end;
        [ct, ax, vc, ind, Nax] = principal_ellipsoid(XXk,center,ax,vec);
        idx = ~ind;
        figure(3);
        plot3(XXk(1,ind),XXk(3,ind),-XXk(2,ind), '.', XXk(1,idx),XXk(3,idx),-XXk(2,idx), '.');
        hold on;
        % generate ellipsoid from unit sphere
        XYZe = vc*diag(ax)*XYZs+ct(:,ones(1,npts));
        Xe(:) = XYZe(1,:);
        Ye(:) = XYZe(2,:);
        Ze(:) = XYZe(3,:);
        surf(Xe,Ze,-Ye,'FaceColor',palette(1,:), 'EdgeColor',palette(1,:), 'FaceAlpha',0.05);
        hold off;
        view(az,el);
        fprintf(1,'\nCheck if the principal ellipsoid need to be resized:\nPress any key to continue...\n');
        figure(3); axis equal off;
        pause;
        re_edit = input('Need to resize the ellipsoid or not? ([]=yes, other=no) ','s');
        re_edit = isempty(re_edit);
    end;
    vec = vc;
    temp = uint8(false(size(ind)));     % idx=~ind are treated as outliners (seg_id == 0)
    temp(ind) = ii;
    seg_id(indk) = temp;
    member_nax(ii,kk) = Nax;
    member_center(:,ii,kk) = ct;
    member_axes(:,ii,kk) = ax;
    
    % project wing points between two planes of a faction of span distance
    % on to the plane perpendicular to span axis and through the center
    XXk = XXk(:,ind);
    Nb = size(XXk,2);
    vc = vc(:,1);
    ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*fracw/Nax;
    idx = vc'*(XXk-ctt(:,1)*ones(1,Nb))<=0 & vc'*(XXk-ctt(:,2)*ones(1,Nb))>=0;
    ind = ~idx;
    figure(3);
    plot3(XXk(1,ind),XXk(3,ind),-XXk(2,ind), '.', XXk(1,idx),XXk(3,idx),-XXk(2,idx), '.');
    hold on;
    XXk = XXk(:,idx);
    Nb = size(XXk,2);
    XXk = XXk-ct(:,ones(1,Nb));
    XXk = XXk-vc*(vc'*XXk);
    if pca_sw,
        [~,~,vc] = svd(XXk*XXk');
        vec(:,2) = vc(:,1);
    else
        [~, vc] = max_distance(XXk);
        vec(:,2) = vc/norm(vc);
    end;
    % plot the two semi axes (span and chord)
    ctt = ct(:,ones(1,2));
    xyz = ctt+vec(:,1:2)*diag(ax(1:2))*nlen_vec;
    for i=1:2,
        plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
    end;
    axis equal vis3d off;
    view(az,el); hold off;
    
    % adjust the orientation of member axes
    fprintf(1, 'Identify the axes orientation (ax3=cross(ax1,ax2)) of member %d:\n',ii);
    fprintf(1,['\nThis function take span (long) direction as x axis and chord direction as y axis!\n' ...
        '\nCheck the direction of axis 1 (blue axis):\nPress any key to continue...\n']);
    pause;
    flag = input('Does the blue axis point to ''+'' or ''-'' direction? ([]=''+'', other=''-'') ','s');
    if ~isempty(flag),
        vec(:,1) = -vec(:,1);
    end;
    fprintf(1,'\nCheck the direction of axis 2 (green axis):\nPress any key to continue...\n');
    pause;
    flag = input('Does the green axis point to ''+'' or ''-'' direction? ([]=''+'', other=''-'') ','s');
    if ~isempty(flag),
        vec(:,2) = -vec(:,2);
    end;
    % rectify axes to right handed
    vec(:,3) = cross(vec(:,1), vec(:,2));
    member_omega(:,ii,kk) = rodrigues(vec);
end;
% save segment index
seg_cell{bk} = seg_id;

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
    ax = member_axes(:,ii,kk);
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
T_step = zeros(3,nparts+1);
Q_step = zeros(4,nparts+1);
Q_step(4,:) = 1;
kk0 = kk;
for kk = Nstart+1:n_frame,
    if kk == base+nfpu+1,
        base = base+nfpu;
        count = count+1;
        nc = sprintf(['%0' ndigit 'd'],count);
        % load visual hull data
        load([imgdir '/3D_points_part' nc '.mat']);
        seg_cell = cell(size(X_cell));
    end;
    if active_images(kk),
        bk = kk-base;
        XX = X_cell{bk};
        
        % initial position
        center0 = member_center(:,:,kk0);
        center = center0+T_step;
        semiax = member_axes(:,:,kk0);
        %     omega = member_omega(:,:,kk0);
        Q0 = trans_quat_axis(member_omega(:,:,kk0));
        Qt = quatmul(Q_step, Q0, 1);
        % segment main body, compute body axes
        ct = center(:,1);
        ax = semiax(:,1);
        %     vec = rodrigues(omega(:,1));
        vec = trans_quat_mat(Qt(:,1));
        [ct, ax, vec, ind, Nax] = principal_ellipsoid(XX,ct,ax,vec,10);
        indk = ~ind;
        seg_id = uint8(ind);
        member_nax(1,kk) = Nax;
        member_center(:,1,kk) = ct;
        member_axes(:,1,kk) = ax;
        
        % refine the heading direction
        vc =  vec(:,1);
        ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*(1-fracb)/Nax;
        XXk = XX(:,ind);
        Nb = size(XXk,2);
        idx = vc'*(XXk-ctt(:,1)*ones(1,Nb))>=0;
        ind = vc'*(XXk-ctt(:,2)*ones(1,Nb))<=0;
        vc = mean(XXk(:,idx),2)-mean(XXk(:,ind),2);
        % Schmidt orthogonalization
        vc = vc/norm(vc);
        vec(:,1) = vc;
        vc = vec(:,2)-vec(:,2)'*vc*vc;
        vec(:,2) = vc/norm(vc);
        vec(:,3) = cross(vec(:,1), vec(:,2));
        member_omega(:,1,kk) = rodrigues(vec);
        
        % compute wing axes
        XXk = XX(:,indk);
        flag = 0;
        for ii=2:nparts+1,
            for jj=ii+1:nparts+1,
                %             intp = ellipsoids_intersection(center(:,[ii,jj]), semiax(:,[ii,jj]), [rodrigues(omega(:,ii)), rodrigues(omega(:,jj))], XXk);
                intp = ellipsoids_intersection(center(:,[ii,jj]), semiax(:,[ii,jj]), [trans_quat_mat(Qt(:,ii)), trans_quat_mat(Qt(:,jj))], XXk);
                if max(intp)>thresh_dist,
                    center(:,2:end) = center0(:,2:end);
                    Qt(:,2:end) = Q0(:,2:end);
                    flag = 1;
                    break;
                end;
            end;
            if flag, break; end;
        end;
        
        temp = uint8(false(1,size(XXk,2)));
        for ii=2:nparts+1,
            ct = center(:,ii);
            ax = semiax(:,ii);
            %         vect = rodrigues(omega(:,ii));
            vect = trans_quat_mat(Qt(:,ii));
            [ct, ax, vec, ind, Nax] = principal_ellipsoid(XXk,ct,ax,vect,10);
            temp(ind) = ii;
            member_nax(ii,kk) = Nax;
            member_center(:,ii,kk) = ct;
            member_axes(:,ii,kk) = ax;
            
            % project wing points between two planes on to the plane perpendicular to span axis
            Xk = XXk(:,ind);
            Nb = size(Xk,2);
            vc = vec(:,1);
            ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*fracw/Nax;
            ind = vc'*(Xk-ctt(:,1)*ones(1,Nb))<=0  & vc'*(Xk-ctt(:,2)*ones(1,Nb))>=0;
            Xk = Xk(:,ind);
            Nb = size(Xk,2);
            Xk = Xk-ct(:,ones(1,Nb));
            Xk = Xk-vc*(vc'*Xk);
            if pca_sw,
                [~,~,vc] = svd(Xk*Xk');
                vc = vc(:,1);
            else
                [~, vc] = max_distance(Xk);
                vc = vc/norm(vc);
            end;
            if vc'*vect(:,2)<0,
                vec(:,2) = -vc;
            else
                vec(:,2) = vc;
            end;
            vec(:,3) = cross(vec(:,1), vec(:,2));
            member_omega(:,ii,kk) = rodrigues(vec);
        end;
        seg_id(indk) = temp; % seg_id == 0 are treated as outliners
        seg_cell{bk} = seg_id;
        
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
            ax = member_axes(:,ii,kk);
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
        
        % update speed
        T_step = member_center(:,:,kk)-member_center(:,:,kk0);
        Q0(1:3,:) = -Q0(1:3,:);
        Q_step = quatmul(trans_quat_axis(member_omega(:,:,kk)), Q0, 1);
        kk0 = kk;
    end;
    if kk == base+nfpu || kk == n_frame,
        fprintf(1,'Saving segmentation result ''segment_part%s.mat''.\n',nc);
        save_name = [imgdir '/segment_part' nc '.mat'];
        eval(['save ' save_name ' X_cell seg_cell']);
    end;
end;
