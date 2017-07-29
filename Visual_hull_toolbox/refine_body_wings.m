kk = track_list(1)-1;
count = ceil(kk/nfpu);
base = (count-1)*nfpu;
nc = sprintf(['%0' ndigit 'd'],count);
save_name = [imgdir '/segment_part' nc '.mat'];
if exist(save_name,'file')==2,
    load(save_name);
else
    load([imgdir '/3D_points_part' nc '.mat']);
end;

kn=0;
for kk = track_list,
    kn = kn+1;
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
        [ct, ax, vec, ind] = ellipsoid_fitting(XXk,center,axes_mean(:,1),vect, MaxIter);
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
            vect = trans_quat_mat(Qt(:,ii-1,kn));
            indk = (seg_id==ii) | (seg_id==0);
            XXk = XX(:, indk);
            [ct, ax, vec, ind] = ellipsoid_rotation2(XXk,root(:,ii-1),axes_mean(:,ii),vect, MaxIter, fracw);
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