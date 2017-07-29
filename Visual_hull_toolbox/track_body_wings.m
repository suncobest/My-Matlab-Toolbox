% track and segment members given initial segmentation
if ~exist('period_sw', 'var') || isempty(period_sw),
    fprintf(1,['\nThis script will apply periodic constraint to:\n'...
        '1: estimate all wing''s direction including axes pointing and incidence;\n'...
        '2: only estimate wing''s incidence while track wing''s heading frame by frame;\n'...
        'Please choose:\n']);
    period_sw = input('Apply periodic constraint to estimate wing''s: ([]=all direction, other=only incidence) ','s');
    period_sw = isempty(period_sw);
end;
if ~exist('func_flag', 'var') || isempty(func_flag),
     flag = 1;
    while flag,
        func_flag = input('Which function to call: (1=ellipsoid_tracking, 2=ellipsoid_rotation, 3=ellipsoid_rotation2) ');
        if length(func_flag)~=1 || all(func_flag~=[1 2 3]),
            fprintf('Unexpected input! Please enter again!\n')
            continue;
        end;
        flag = 0;
    end;
end;

kk = track_list(1)-1;
count = ceil(kk/nfpu);
base = (count-1)*nfpu;
nc = sprintf(['%0' ndigit 'd'],count);
load([imgdir '/3D_points_part' nc '.mat']);

kk0 = find(active_images(1:kk),1,'last');
for kk = track_list,
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
        center = member_center(:,1,kk0);
        omega = member_omega(:,:,kk0);
        R = rodrigues(omega(:,1));
        % segment main body, compute body axes
        [ct, ax, vec, ind] = ellipsoid_fitting(XX,center,axes_mean(:,1),R,MaxIter);
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
        % Schmidt orthogonalization
        body_x = vc/norm(vc);
        
        % compute wing axes
        root = R*wroot+ct(:,[1 1]);
        for ii=2:nparts+1,
            indk = ~seg_id;    % segment points 0
            XXk = XX(:,indk);
            temp = uint8(false(1,size(XXk,2)));
            if period_sw,
                % apply periodic constraint to estimate wing's direction (including axis and incidence)
                vect = R*trans_quat_mat(Qbw(:,ii-1,kk));
            else
                % track wing's axis and apply periodic constraint to estimate wing's incidence
                vect = rodrigues(omega(:,ii));
                vc = R*trans_quat_mat(Qbw(:,ii-1,kk));
                vc = vc(:,2)-vc(:,2)'*vect(:,1)*vect(:,1);
                vect(:,2) = vc/norm(vc);
                vect(:,3) = cross(vect(:,1),vect(:,2));
            end
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
        kk0 = kk;
    end;
    if kk == base+nfpu || kk == n_frame,
        Nb = sum(active_images(1:kk));
        wroot = sum(root_mat(:,:,active_images(1:kk)),3)/Nb;
        wroot(2,2) = -wroot(2,2);
        wroot = sum(wroot,2)/2*[1,1];
        wroot(2,2) = -wroot(2,2);
        fprintf(1,'Saving segmentation result ''segment_part%s.mat''.\n',nc);
        save_name = [imgdir '/segment_part' nc '.mat'];
        eval(['save ' save_name ' X_cell seg_cell']);
    end;
end;