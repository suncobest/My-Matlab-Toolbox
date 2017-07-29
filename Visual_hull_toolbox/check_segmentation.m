check_direction = input('Need to check member direction of every frame? ([]=yes, other=no) ','s');
check_direction = isempty(check_direction);

Nb = frame_list(1);
kk = track_list(1);
count = ceil(kk/nfpu);
base = (count-1)*nfpu;
nc = sprintf(['%0' ndigit 'd'],count);
save_name = [imgdir '/segment_part' nc '.mat'];
if exist(save_name,'file')==2,
    load(save_name);
else
    if count>1,
        fprintf(1,'\nERROR: cannot find data ''%s''!\n',save_name);
        return;
    else
        assert(frame_list(end)<nfpu, 'All frames are assumed in part 1!');
        load([imgdir '/3D_points_part' nc '.mat']);
    end;
end;

for kk = track_list,
    if kk == base+nfpu+1,
        base = base+nfpu;
        count = count+1;
        nc = sprintf(['%0' ndigit 'd'],count);
        % load visual hull data
        save_name = [imgdir '/segment_part' nc '.mat'];
        if exist(save_name,'file')==2,
            load(save_name);
        else
            fprintf(1,'\nERROR: cannot find data ''%s''!\n',save_name);
            return;
        end;
    end;
    if active_images(kk),
        bk = kk-base;
        XX = X_cell{bk};
        seg_id = seg_cell{bk};
        figure(3); hold off;
        idx = (seg_id==0);
        ii = nparts+2;
        plot3(XX(1,idx),XX(3,idx),-XX(2,idx), '.', 'color',palette(ii,:));
        hold on;
        for ii=1:nparts+1,
            idx = (seg_id==ii);
            plot3(XX(1,idx),XX(3,idx),-XX(2,idx), '.', 'color',palette(ii,:));
            center = member_center(:,ii,kk);
            semiax = member_axes(:,ii,kk);
            vec = rodrigues(member_omega(:,ii,kk));
            ctt = center(:,ones(1,2));
            xyz = ctt+vec(:,1:2)*diag(semiax(1:2))*nlen_vec;
            for i=1:2,
                plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
            end;
        end;
        title(['Segment result of frame: ', num2str(kk)]);
        axis equal vis3d off;
        view(az,el); hold off;
        pause(0.01);
        if check_direction,
            fprintf(1,'\nCheck axes direction of frame %d:\n',kk);
            flag = input('Some parts'' direciton need to be rectified? ([]=no, other=yes) ','s');
            if isempty(flag),
                Not_ok(kk-Nb+1) = 0;
            end;
        end;
    end;
end;

if ~check_direction,
    flag = input('Segmentation of all frames look all right? ([]=yes, other=no)');
    Not_ok(:) = 0;
    if ~isempty(flag),
        ind = input(['Enter number of frames (' num2str(track_list(1)) '~' num2str(track_list(end)) ') that need to be rectified: ']);
        Not_ok(ind-Nb+1) = 1;
    end;
end;