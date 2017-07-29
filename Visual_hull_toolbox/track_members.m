% track and segment members given initial segmentation
kk = track_list(1);
if ~active_images(kk),
    fprintf(1,'\nProgram need the start frame to be active!\n');
    return;
end;
count = ceil(kk/nfpu);
base = (count-1)*nfpu;
bk = kk-base;
nc = sprintf(['%0' ndigit 'd'],count);
load([imgdir '/3D_points_part' nc '.mat']);
XX = X_cell{bk};

if kk==1,
    seg_id = seg_cell{bk};
    root = NaN(3,nparts);

    MaxIter = input('Maximum iteration number for the ellipsoid_tracking function: ([]=20) ');
    if isempty(MaxIter),
        MaxIter = 20;
    end;
    fracw = input('Fraction number for wings to estimate chord direction: ([]=0.1) ');
    if isempty(fracw),
        fracw = 0.1;
    end;
    fracb = input('Fraction number for body to estimate heading direction: ([]=0.2) ');
    if isempty(fracb),
        fracb = 0.2;
    end;
    
    % compute body axes
    center = member_center(:,1,kk);
    axes_mean = member_axes(:,:,kk);
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
            semiax = axes_mean(:,1);
            fprintf(1,'Do not resize at all!\n');
        elseif (n~=1 && n~=3) || any(Nb<=0),
            fprintf(1,'Unexpected input! Please enter again!\n');
            continue;
        else
            semiax = axes_mean(:,1).*Nb(:);
        end;
        [ct, ax, vc, ind] = ellipsoid_fitting(XXk,center,semiax,vec);
        idx = ~ind;
        figure(3);
        plot3(XXk(1,ind),XXk(3,ind),-XXk(2,ind), '.', XXk(1,idx),XXk(3,idx),-XXk(2,idx), '.');
        hold on;
        % generate ellipsoid from unit sphere
        XYZe = vc*diag(semiax)*XYZs+ct(:,ones(1,npts));
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
    member_center(:,1,kk) = ct;
    member_axes(:,1,kk) = ax;
    axes_mean(:,1) = semiax;
    
    % refine the heading direction
    vc =  vec(:,1);
    ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*(1-fracb);
    ax0 = XXk(:,ind);
    Nb = size(ax0,2);
    idx = vc'*(ax0-ctt(:,1)*ones(1,Nb))>=0;
    ind = vc'*(ax0-ctt(:,2)*ones(1,Nb))<=0;
    vc = mean(ax0(:,idx),2)-mean(ax0(:,ind),2);
    % Schmidt orthogonalization
    vc = vc/norm(vc);
    vec(:,1) = vc;
    body_x = vc;
    vc = vec(:,2)-vec(:,2)'*vc*vc;
    vec(:,2) = vc/norm(vc);
    vec(:,3) = cross(vec(:,1), vec(:,2));
    
    indk = ~(ind | idx);
    figure(3);
    plot3(ax0(1,indk),ax0(3,indk),-ax0(2,indk), '.', ax0(1,ind),ax0(3,ind),-ax0(2,ind), '.', ax0(1,idx),ax0(3,idx),-ax0(2,idx), '.');
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
        fprintf(1,'\nGenerate segment ellipsoid for part %d (%s):\n',ii,color_cell{ii});
        % generate ellipsoid from unit sphere
        indk = (seg_id==ii);
        XXk = XX(:, indk);
        [center, semiax0, vect] = equivalent_ellipsoid(XXk);
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
                semiax = semiax0;
                fprintf(1,'Do not resize at all!\n');
            elseif (n~=1 && n~=3) || any(Nb<=0),
                fprintf(1,'Unexpected input! Please enter again!\n');
                continue;
            else
                semiax = semiax0.*Nb(:);
            end;
            [ct, ax, vec, ind] = ellipsoid_fitting(XXk,center,semiax,vect);
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
            fprintf(1,'\nCheck if the principal ellipsoid need to be resized:\n');
            figure(3); axis equal off;
            re_edit = input('Need to resize the ellipsoid or not? ([]=yes, other=no) ','s');
            re_edit = isempty(re_edit);
        end;
        temp = uint8(false(size(ind)));     % idx=~ind are treated as outliners (seg_id == 0)
        temp(ind) = ii;
        seg_id(indk) = temp;
        member_center(:,ii,kk) = ct;
        member_axes(:,ii,kk) = ax;
        axes_mean(:,ii) = semiax;
        
        % project wing points between two planes of a faction of span distance
        % on to the plane perpendicular to span axis and through the center
        XXk = XXk(:,ind);
        Nb = size(XXk,2);
        vc = vec(:,1);
        ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*fracw;
        idx = vc'*(XXk-ctt(:,1)*ones(1,Nb))<=0 & vc'*(XXk-ctt(:,2)*ones(1,Nb))>=0;
        ind = ~idx;
        Xk = XXk(:,idx);
        Xk = Xk-ct(:,ones(1,size(Xk,2)));
        Xk = Xk-vc*(vc'*Xk);
        % choose method to measure chord direction
        pca_sw = input('Which method to estimate wings'' chord direction: ([]=PCA, other=Max distance) ','s');
        pca_sw = isempty(pca_sw);
        re_edit = 1;
        while re_edit,
            figure(3);
            plot3(ax0(1,:),ax0(3,:),-ax0(2,:), '.', XXk(1,ind),XXk(3,ind),-XXk(2,ind), '.', XXk(1,idx),XXk(3,idx),-XXk(2,idx), '.');
            hold on;
            if pca_sw,
                title('PCA');
                [~,~,vc] = svd(Xk*Xk');
                vec(:,2) = vc(:,1);
            else
                title('Max distance');
                [~, vc] = max_projection_vector(Xk);
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
            re_edit = input('Reset method or not? ([]=yes, other=no) ','s');
            re_edit = isempty(re_edit);
            if re_edit,
                pca_sw = ~pca_sw;
            end;
        end;
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
        % compute wing root
        root(:,ii-1) = ct-vec(:,1)*ax(1);
    end;
    axes_mean(:,2:3) = sum(axes_mean(:,2:3),2)/2*[1,1];
    
else
    
    kk0 = find(active_images(1:kk-1),1,'last');
    % initial position
    center = member_center(:,1,kk0);
    omega = member_omega(:,:,kk0);
    % segment main body, compute body axes
    vect = rodrigues(omega(:,1));
    [ct, ax, vec, ind] = ellipsoid_fitting(XX,center,axes_mean(:,1),vect,MaxIter);
    seg_id = uint8(ind);
    member_center(:,1,kk) = ct;
    member_axes(:,1,kk) = ax;
    
    % refine the heading direction
    vc =  vec(:,1);
    ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*(1-fracb);
    ax0 = XX(:,ind);
    Nb = size(ax0,2);
    idx = vc'*(ax0-ctt(:,1)*ones(1,Nb))>=0;
    ind = vc'*(ax0-ctt(:,2)*ones(1,Nb))<=0;
    vc = mean(ax0(:,idx),2)-mean(ax0(:,ind),2);
    body_x = vc/norm(vc);
    
    % compute wing's roots
    root = vect*wroot+ct(:,[1 1]);
    for ii=2:nparts+1,
        indk = ~seg_id;    % segment points 0
        XXk = XX(:,indk);
        temp = uint8(false(1,size(XXk,2)));
        vect =  rodrigues(omega(:,ii));
        [ct, ax, vec, ind] = ellipsoid_rotation(XXk,root(:,ii-1),axes_mean(:,ii),vect,MaxIter);
        temp(ind) = ii;
        seg_id(indk) = temp;    % seg_id == 0 are treated as outliners
        member_center(:,ii,kk) = ct;
        member_axes(:,ii,kk) = ax;
        
        % project wing points between two planes of a faction of span distance
        % on to the plane perpendicular to span axis and through the center
        XXk = XXk(:,ind);
        Nb = size(XXk,2);
        vc = vec(:,1);
        ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*fracw;
        idx = vc'*(XXk-ctt(:,1)*ones(1,Nb))<=0 & vc'*(XXk-ctt(:,2)*ones(1,Nb))>=0;
        ind = ~idx;
        Xk = XXk(:,idx);
        Xk = Xk-ct(:,ones(1,size(Xk,2)));
        Xk = Xk-vc*(vc'*Xk);
        vect = vect(:,2)-vc'*vect(:,2)*vc;
        % choose method to measure chord direction
        re_edit = 1;
        while re_edit,
            flag = 1;
            while flag,
                fprintf(1,'\nMethods: []=1=PCA; 2=Max distance; 3=Max projection;\n');
                pca_sw = input('Which method to estimate wings'' chord direction:');
                if isempty(pca_sw),
                    pca_sw = 1;
                end;
                switch pca_sw,
                    case 1,
                        title('PCA');
                        [~,~,vc] = svd(Xk*Xk');
                        vec(:,2) = vc(:,1);
                    case 2,
                        title('Max distance');
                        [~, vc] = max_projection_vector(Xk);
                        vec(:,2) = vc/norm(vc);
                    case 3,
                        title('Max projection');
                        [~, vc] = max_projection_vector(Xk,vect);
                        vec(:,2) = vc/norm(vc);
                    otherwise,
                        disp('Uexpected input! Please enter again!')
                        continue;
                end;
                flag = 0;
            end;
            figure(3);
            plot3(ax0(1,:),ax0(3,:),-ax0(2,:), '.', XXk(1,ind),XXk(3,ind),-XXk(2,ind), '.', XXk(1,idx),XXk(3,idx),-XXk(2,idx), '.');
            hold on;
            % plot the two semi axes (span and chord)
            ctt = ct(:,ones(1,2));
            xyz = ctt+vec(:,1:2)*diag(ax(1:2))*nlen_vec;
            for i=1:2,
                plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
            end;
            axis equal vis3d off;
            view(az,el); hold off;
            re_edit = input('Reset method or not? ([]=yes, other=no) ','s');
            re_edit = isempty(re_edit);
        end;
        
        % adjust the orientation of member axes
        fprintf(1, 'Identify the axes orientation (ax3=cross(ax1,ax2)) of member %d:\n',ii);
        fprintf(1,'\nCheck the direction of axis 2 (green axis):\n');
        flag = input('Does the green axis point to ''+'' or ''-'' direction? ([]=''+'', other=''-'') ','s');
        if ~isempty(flag),
            vec(:,2) = -vec(:,2);
        end;
        % rectify axes to right handed
        vec(:,3) = cross(vec(:,1), vec(:,2));
        member_omega(:,ii,kk) = rodrigues(vec);
        % compute wing root
        root(:,ii-1) = ct-vec(:,1)*ax(1);
    end;
end;

% save segment index
seg_cell{bk} = seg_id;

% compute lateral direction of body (from left to right wing center)
% body_y = root(:,1)-root(:,2);
body_y = member_center(:,2,kk)-member_center(:,3,kk);
body_y = body_y-body_y'*body_x*body_x;     % orthogonalization
body_y = body_y/norm(body_y);   %  axis2 of body
% update the body orientation
vec = [body_x, body_y, cross(body_x,body_y)];
member_omega(:,1,kk) = rodrigues(vec);
root_mat(:,:,kk) = vec'*(root-member_center(:,1,kk)*[1 1]);
Nb = sum(active_images(1:kk));
wroot = sum(root_mat(:,:,active_images(1:kk)),3)/Nb;
wroot(2,2) = -wroot(2,2);
wroot = sum(wroot,2)/2*[1,1];
wroot(2,2) = -wroot(2,2);

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
if kk == base+nfpu || kk == n_frame,
    fprintf(1,'Saving segmentation result ''segment_part%s.mat''.\n',nc);
    save_name = [imgdir '/segment_part' nc '.mat'];
    eval(['save ' save_name ' X_cell seg_cell']);
end;

if kk==track_list(end),
    return;
end;

for kk = track_list(2:end),
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
        kk0 = find(active_images(1:kk-1),1,'last');
        center = member_center(:,1,kk0);
        omega = member_omega(:,:,kk0);
        % segment main body, compute body axes
        vect = rodrigues(omega(:,1));
        [ct, ax, vec, ind] = ellipsoid_fitting(XX,center,axes_mean(:,1),vect,MaxIter);
        seg_id = uint8(ind);
        member_center(:,1,kk) = ct;
        member_axes(:,1,kk) = ax;
        
        % refine the heading direction
        vc =  vec(:,1);
        ctt = ct(:,[1,1])+vc*[1,-1]*ax(1)*(1-fracb);
        XXk = XX(:,ind);
        Nb = size(XXk,2);
        idx = vc'*(XXk-ctt(:,1)*ones(1,Nb))>=0;
        ind = vc'*(XXk-ctt(:,2)*ones(1,Nb))<=0;
        vc = mean(XXk(:,idx),2)-mean(XXk(:,ind),2);
        body_x = vc/norm(vc);
        
        % compute wing's roots
        root = vect*wroot+ct(:,[1 1]);
        for ii=2:nparts+1,
            indk = ~seg_id;    % segment points 0
            XXk = XX(:,indk);
            temp = uint8(false(1,size(XXk,2)));
            vect = rodrigues(omega(:,ii));
            [ct, ax, vec, ind] = ellipsoid_tracking(XXk,root(:,ii-1),axes_mean(:,ii),vect,MaxIter, fracw);
            temp(ind) = ii;
            seg_id(indk) = temp;    % seg_id == 0 are treated as outliners
            member_center(:,ii,kk) = ct;
            member_axes(:,ii,kk) = ax;
            member_omega(:,ii,kk) = rodrigues(vec);
            % compute wing root
            root(:,ii-1) = ct-vec(:,1)*ax(1);
        end;
        seg_cell{bk} = seg_id;
        % compute lateral direction of body (from left to right wing joint)
        % body_y = root(:,1)-root(:,2);
        body_y = member_center(:,2,kk)-member_center(:,3,kk);
        body_y = body_y-body_y'*body_x*body_x;     % orthogonalization
        body_y = body_y/norm(body_y);   %  axis2 of body
        % update the body orientation
        vec = [body_x, body_y, cross(body_x,body_y)];
        member_omega(:,1,kk) = rodrigues(vec);
        root_mat(:,:,kk) = vec'*(root-member_center(:,1,kk)*[1 1]);
        Nb = sum(active_images(1:kk));
        wroot = sum(root_mat(:,:,active_images(1:kk)),3)/Nb;
        wroot(2,2) = -wroot(2,2);
        wroot = sum(wroot,2)/2*[1,1];
        wroot(2,2) = -wroot(2,2);
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

if kk>nfpu && mod(kk,nfpu) && kk~=n_frame,
    fprintf(1,'Saving segmentation result ''segment_part%s.mat''.\n',nc);
    save_name = [imgdir '/segment_part' nc '.mat'];
    eval(['save ' save_name ' X_cell seg_cell']);
end;
