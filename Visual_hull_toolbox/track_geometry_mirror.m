% track_geometry_mirror
% regionprops uses 8-connectivity to compute binary images, (calls bwconncomp inside function)
% so I need to call bwconncomp first to use 4-connectivity to compute Area,
% BoundingBox, and Centroid of each region

frame_kk = imread([imPreNumf strnum_frame{kk} '.png']);
geom_kk = regionprops(bwconncomp(frame_kk,4),'basic');
% geometric feature
boundbox_kk = vertcat(geom_kk.BoundingBox)';
centroid_kk = vertcat(geom_kk.Centroid)';
area_kk = horzcat(geom_kk.Area);

threshold_area = max(area_kk, [], 2)*ndarea;
nblob = length(area_kk);

figure(2);
image(frame_kk);
title(['Check camera number in every view of frame ' num2str(kk) ':']);
colormap(gray(2));      % colormap for binary image
axis image;
hold on;
active_view = active_imgviews(:,kk);
for pp = 1:n_cam,
    if active_view(pp),
        kth = (kk-1)*n_cam+pp;
        xx = centroid_last(:,pp);
        % find the closest center of blobs to the last center postion
        [xx, ind] = min(sum((centroid_kk-repmat(xx,1,nblob)).^2, 1), [], 2);
        center = centroid_kk(:,ind);
        if sqrt(xx)>distance_toler(pp) || area_kk(ind) < threshold_area,
            fprintf(1,'\nThe foreground of (camera %d, frame %d) is abnormal!\n',pp,kk);
            active_view(pp) = 0;
        else
            temp = boundbox_kk(:,ind);
            bounding_mat(:,kth) = temp;
            centroid_last(:,pp) = center;
            center2d_mat(:,kth) = center;
            area_mat(:,kth) = area_kk(ind);
            temp(1:2) = ceil(temp([2 1]));
            foreground_cell{kth} = frame_kk(temp(1) : temp(1)+temp(4)-1, temp(2) : temp(2)+temp(3)-1);
            % plot camera number
            text(center(1),center(2),['\it\fontname{Arial}\fontsize{9}\color{green}' ...
                        num2str(pp)],'HorizontalAlignment','center');
            rectangle('Position',boundbox_kk(:,ind),'edgecolor','b','linewidth',1);
        end;
    end;
end;
hold off;
active_imgviews(:,kk) = active_view;