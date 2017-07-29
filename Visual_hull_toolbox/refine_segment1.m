% track and segment members given initial segmentation
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
save_name = [imgdir '/segment_part' nc '.mat'];
if exist(save_name,'file')==2,
    load(save_name);
else
    fprintf(1,'\nERROR: cannot find data ''%s''!\n',save_name);
    return;
end;
XX = X_cell{bk};
seg_id = seg_cell{bk};

MaxIter = input('Maximum iteration number for the ellipsoid_fitting function: ([]=50) ');
if isempty(MaxIter),
    MaxIter = 50;
end;

% Default: do not optimize orientation of ellipsoid
SW = input('Optimize orientation or not? ([]=no (recommended), other=yes) ', 's');
SW = ~isempty(SW);

% compute body axes
center = member_center(:,1,kk);
semiax0 = axes_mean(:,1);
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
        fprintf(1,'Do not resize at all!\n');
        semiax = semiax0;
    elseif (n~=1 && n~=3) || any(Nb<=0),
        fprintf(1,'Unexpected input! Please enter again!\n');
        continue;
    else
        semiax = semiax0.*Nb(:);
    end;
    [ct, ax, vc, ind] = ellipsoid_fitting(XXk,center,semiax,vec, MaxIter, SW);
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
vc = vec(:,2)-vec(:,2)'*vc*vc;
vec(:,2) = vc/norm(vc);
vec(:,3) = cross(vec(:,1), vec(:,2));
member_omega(:,1,kk) = rodrigues(vec);

% compute wing axes
for ii=2:nparts+1,
    center = member_center(:,ii,kk);
    semiax0 = axes_mean(:,ii);
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
            semiax = semiax0;
        elseif (n~=1 && n~=3) || any(Nb<=0),
            fprintf(1,'Unexpected input! Please enter again!\n');
            continue;
        else
            semiax = semiax0.*Nb(:);
        end;
        [ct, ax, vec, ind] = ellipsoid_fitting(XXk,center,semiax,vect, MaxIter, SW);
        idx = ~ind;
        figure(3);
        plot3(XXk(1,ind),XXk(3,ind),-XXk(2,ind), '.', XXk(1,idx),XXk(3,idx),-XXk(2,idx), '.');
        hold on;
        % generate ellipsoid from unit sphere
        XYZe = vec*diag(ax)*XYZs+ct(:,ones(1,npts));
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
    temp = uint8(false(size(ind)));     % idx=~ind are treated as outliners (seg_id == 0)
    temp(ind) = ii;
    seg_id(indk) = temp;
    member_center(:,ii,kk) = ct;
    member_axes(:,ii,kk) = ax;
    axes_mean(:,ii) = semiax;    % update segment ellipsoid
    
    % project wing points between two planes of a faction of span distance
    % on to the plane perpendicular to span axis and through the center
    XXk = XXk(:,ind);
    Nb = size(XXk,2);
    vc = vec(:,1);
    ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*fracw;
    idx = vc'*(XXk-ctt(:,1)*ones(1,Nb))<=0 & vc'*(XXk-ctt(:,2)*ones(1,Nb))>=0;
    XXk = XXk(:,idx);
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
        semiax = axes_mean(:,1);
        vec = rodrigues(member_omega(:,1,kk));
        indk = (seg_id==1) | (seg_id==0);
        XXk = XX(:, indk);
        [ct, ax, vec, ind] = ellipsoid_fitting(XXk,center,semiax,vec, MaxIter, SW);
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
        % Schmidt orthogonalization
        vc = vc/norm(vc);
        vec(:,1) = vc;
        vc = vec(:,2)-vec(:,2)'*vc*vc;
        vec(:,2) = vc/norm(vc);
        vec(:,3) = cross(vec(:,1), vec(:,2));
        member_omega(:,1,kk) = rodrigues(vec);
       
        % compute wing axes
        for ii=2:nparts+1,
            center = member_center(:,ii,kk);
            semiax = axes_mean(:,ii);
            vect = rodrigues(member_omega(:,ii,kk));
            indk = (seg_id==ii) | (seg_id==0);
            XXk = XX(:, indk);
            [ct, ax, vec, ind] = ellipsoid_fitting(XXk,center,semiax,vect, MaxIter, SW);
            temp = uint8(false(size(ind)));     % idx=~ind are treated as outliners (seg_id == 0)
            temp(ind) = ii;
            seg_id(indk) = temp;
            member_center(:,ii,kk) = ct;
            member_axes(:,ii,kk) = ax;
            
            % project wing points between two planes on to the plane perpendicular to span axis
            XXk = XXk(:,ind);
            Nb = size(XXk,2);
            vc = vec(:,1);
            ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*fracw;
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
    end;
    if kk == base+nfpu || kk == n_frame,
        fprintf(1,'Saving segmentation result ''segment_part%s.mat''.\n',nc);
        save_name = [imgdir '/segment_part' nc '.mat'];
        eval(['save ' save_name ' X_cell seg_cell']);
    end;
end;
