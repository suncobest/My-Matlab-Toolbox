% track and segment members given initial segmentation
% Nstart = 1;
kk = Nstart;
if ~active_images(kk),
    fprintf(1,'\nProgram need the start frame to be active!\n');
    return;
end;
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

fracw = input('Fraction number for wings to estimate chord direction: ([]=0.05) ');
if isempty(fracw),
    fracw = 0.05;
end;
fracb = input('Fraction number for body to estimate heading direction: ([]=0.2) ');
if isempty(fracb),
    fracb = 0.2;
end;
pca_sw = input('Method to estimate wings'' chord direction: ([]=PCA, other=max_projection_vector) ','s');
pca_sw = isempty(pca_sw);

thresh_dist = input('Distance threshold among members beside body: ([]=0.1) ');
if isempty(thresh_dist),
    thresh_dist = 0.1;
end;

len_chain = input('How many past frames do you want to track? ([]=4) ');
if isempty(len_chain),
    len_chain = 4;
else
    len_chain = round(len_chain);
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
    fprintf(1,'Check if the principal ellipsoid need to be resized:\n');
    figure(3); axis equal off;
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
ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*(1-fracb)/Nax;
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
fprintf(1,'\nCheck the definition of heading direction (blue axis):\n');
pause(2);

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
    fprintf(1,'\nGenerate principal ellipsoid for part %d (%s):\n',ii,color_cell{ii});
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
        fprintf(1,'\nResize section ellipsoid of part %d...\n',ii);
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
        fprintf(1,'\nCheck if the principal ellipsoid need to be resized:\n');
        figure(3); axis equal off;
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
        [~, vc] = max_projection_vector(XXk);
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
        '\nCheck the direction of axis 1 (blue axis):\n']);
    flag = input('Does the blue axis point to ''+'' or ''-'' direction? ([]=''+'', other=''-'') ','s');
    if ~isempty(flag),
        vec(:,1) = -vec(:,1);
    end;
    fprintf(1,'\nCheck the direction of axis 2 (green axis):\n');
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
if kk == base+nfpu || kk == n_frame,
    fprintf(1,'Saving segmentation result ''segment_part%s.mat''.\n',nc);
    save_name = [imgdir '/segment_part' nc '.mat'];
    eval(['save ' save_name ' X_cell seg_cell']);
end;

if Nstart==n_frame,
    return;
end;
Nstop = Nstart+1;
if ~active_images(Nstop),
    fprintf(1,'\nProgram need the frame after start frame to be active!\n');
    return;
end;
if len_chain == 1,
    Nstop = n_frame;
end;

kk0 = kk;
for kk = Nstart+1:Nstop,
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
        center = member_center(:,:,kk0);
        semiax = member_axes(:,:,kk0);
        omega = member_omega(:,:,kk0);
        % segment main body, compute body axes
        ct = center(:,1);
        ax = semiax(:,1);
        vec = rodrigues(omega(:,1));
        [ct, ax, vec, ind, Nax] = principal_ellipsoid(XX,ct,ax,vec,MaxIter);
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
        for ii=2:nparts+1,
            indk = ~seg_id;    % segment points 0
            XXk = XX(:,indk);
            temp = uint8(false(1,size(XXk,2)));
            ct = center(:,ii);
            ax = semiax(:,ii);
            vect = rodrigues(omega(:,ii));
            [ct, ax, vec, ind, Nax] = principal_ellipsoid(XXk,ct,ax,vect,MaxIter);
            temp(ind) = ii;
            seg_id(indk) = temp;    % seg_id == 0 are treated as outliners
            member_nax(ii,kk) = Nax;
            member_center(:,ii,kk) = ct;
            member_axes(:,ii,kk) = ax;
            
            % project wing points between two planes on to the plane perpendicular to span axis
            XXk = XXk(:,ind);
            Nb = size(XXk,2);
            vc = vec(:,1);
            ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*fracw/Nax;
            ind = vc'*(XXk-ctt(:,1)*ones(1,Nb))<=0 & vc'*(XXk-ctt(:,2)*ones(1,Nb))>=0;
            XXk = XXk(:,ind);
            Nb = size(XXk,2);
            XXk = XXk-ct(:,ones(1,Nb));
            XXk = XXk-vc*(vc'*XXk);
            vect = vect(:,2)-vc'*vect(:,2)*vc;
            if pca_sw,
                [~,~,vc] = svd(XXk*XXk');
                vc = vc(:,1);
                if vc'*vect<0,
                    vc = -vc;
                end;
            else
                [~, vc] = max_projection_vector(XXk,vect);
                vc = vc/norm(vc);
            end;
            vec(:,2) = vc;
            vec(:,3) = cross(vec(:,1), vec(:,2));
            member_omega(:,ii,kk) = rodrigues(vec);
        end;
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
        kk0 = kk;
    end;
    if kk == base+nfpu || kk == n_frame,
        fprintf(1,'Saving segmentation result ''segment_part%s.mat''.\n',nc);
        save_name = [imgdir '/segment_part' nc '.mat'];
        eval(['save ' save_name ' X_cell seg_cell']);
    end;
end;

if Nstop == n_frame,
    return;
end;
nn = 3*(nparts+1);
id_chain = [Nstart, Nstop];
% spline and squad interpolation of translation and quaternion rotaion
for kk = Nstop+1:n_frame,
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
        kk0 = id_chain(end);
        center0 = member_center(:,:,kk0);
        Q0 = trans_quat_axis(member_omega(:,:,kk0));
        Qt = Q0;
        semiax = member_axes(:,:,kk0);
        center = reshape(spline_interp3(id_chain, reshape(member_center(:,:,id_chain), nn,[]), kk), 3,[]);
        for ii=1:nparts+1,
            Qt(:,ii) = squad(id_chain, trans_quat_axis(reshape(member_omega(:,ii,id_chain),3,[])), kk);
        end;
        % segment main body, compute body axes
        ct = center(:,1);
        ax = semiax(:,1);
        vec = trans_quat_mat(Qt(:,1));
        [ct, ax, vec, ind, Nax] = principal_ellipsoid(XX,ct,ax,vec,MaxIter);
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
        
        % check if predicted members except body have intersection or not
        XXk = XX(:,~seg_id);
        flag = 0;
        for ii=2:nparts+1,
            for jj=ii+1:nparts+1,
                temp = ellipsoids_intersection(center(:,[ii,jj]), semiax(:,[ii,jj]), [trans_quat_mat(Qt(:,ii)), trans_quat_mat(Qt(:,jj))], XXk);
                if max(temp)>thresh_dist,
                    % use position and orientation of last frame
                    center(:,2:end) = center0(:,2:end);
                    Qt(:,2:end) = Q0(:,2:end);
                    flag = 1;
                    break;
                end;
            end;
            if flag, break; end;
        end;
        
        % compute wing axes
        for ii=2:nparts+1,
            indk = ~seg_id;    % segment points 0
            XXk = XX(:,indk);
            Nb = size(XXk,2);
            temp = uint8(false(1,size(XXk,2)));
            ct = center(:,ii);
            ax = semiax(:,ii);
            vect = trans_quat_mat(Qt(:,ii));
            [ct, ax, vec, ind, Nax] = principal_ellipsoid(XXk,ct,ax,vect,MaxIter);
            % use position and orientation of last frame for bad prediction
            if sum(ind)<Nb*thresh_dist,
                ct = center0(:,ii);
                ax = semiax(:,ii);
                vect = trans_quat_mat(Q0(:,ii));
                [ct, ax, vec, ind, Nax] = principal_ellipsoid(XXk,ct,ax,vect, MaxIter);
            end;
            
            temp(ind) = ii;
            seg_id(indk) = temp; % seg_id == 0 are treated as outliners
            member_nax(ii,kk) = Nax;
            member_center(:,ii,kk) = ct;
            member_axes(:,ii,kk) = ax;
            
            % project wing points between two planes on to the plane perpendicular to span axis
            XXk = XXk(:,ind);
            Nb = size(XXk,2);
            vc = vec(:,1);
            ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*fracw/Nax;
            ind = vc'*(XXk-ctt(:,1)*ones(1,Nb))<=0 & vc'*(XXk-ctt(:,2)*ones(1,Nb))>=0;
            XXk = XXk(:,ind);
            Nb = size(XXk,2);
            XXk = XXk-ct(:,ones(1,Nb));
            XXk = XXk-vc*(vc'*XXk);
            vect = vect(:,2)-vc'*vect(:,2)*vc;
            if pca_sw,
                [~,~,vc] = svd(XXk*XXk');
                vc = vc(:,1);
                if vc'*vect<0,
                    vc = -vc;
                end;
            else
                [~, vc] = max_projection_vector(XXk,vect);
                vc = vc/norm(vc);
            end;
            vec(:,2) = vc;
            vec(:,3) = cross(vec(:,1), vec(:,2));
            member_omega(:,ii,kk) = rodrigues(vec);
        end;
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
        id_chain = [id_chain, kk];
        if length(id_chain)>len_chain,
            id_chain(1) = [];
        end;
    end;
    if kk == base+nfpu || kk == n_frame,
        fprintf(1,'Saving segmentation result ''segment_part%s.mat''.\n',nc);
        save_name = [imgdir '/segment_part' nc '.mat'];
        eval(['save ' save_name ' X_cell seg_cell']);
    end;
end;
