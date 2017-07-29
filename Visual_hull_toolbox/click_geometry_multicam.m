% click_geometry_multicam
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
colormap(gray(2));      % colormap for binary image
axis image;
hold on;
% click one region in case that there are multiple valid blobs
click_view = 1;
kth = (kk-1)*n_cam+pp;
while click_view,
    fprintf(1,'\nPlease click the foreground center of (camera %d, frame %d) ...\n',pp,kk);
    title(['Click the foreground center of (camera ' num2str(pp) ', frame ' num2str(kk) '):']);
    plot(centroid_kk(1,:),centroid_kk(2,:),'r+');
    xx = ginput(1);
    % find the closest center of blobs to the clicked postion
    [~, ind] = min(sum((centroid_kk-repmat(xx',1,nblob)).^2, 1), [], 2);
    center = centroid_kk(:,ind);
    % Just click where xx is 0 or out of image if the view is not available
    xx = round(xx);
    flag = xx(1)<0.5 | xx(1)>nx+0.5 | xx(2)<0.5 | xx(2)>ny+0.5;
    if flag || ~frame_kk(xx(2),xx(1)) || area_kk(ind) < threshold_area,
        fprintf(1,'The foreground of (camera %d, frame %d) is abnormal!\n',pp,kk);
        click_view = input('Is this view valid or not?? ([]=no, other=yes) ','s');
        click_view = ~isempty(click_view);
        if ~click_view,
            fprintf(1,'(camera %d, frame %d) is set inactive!\n',pp,kk);
            active_view(kk) = 0;
        end;
    else
        click_view = 0;
        temp = boundbox_kk(:,ind);
        bounding_mat(:,kth) = temp;
        center2d_mat(:,kth) = center;
        area_mat(:,kth) = area_kk(ind);
        temp(1:2) = ceil(temp([2 1]));
        foreground_cell{kth} = frame_kk(temp(1) : temp(1)+temp(4)-1, temp(2) : temp(2)+temp(3)-1);
        % plot camera number
        text(center(1),center(2),['\it\fontname{Arial}\fontsize{9}\color{green}' ...
                        num2str(pp)],'HorizontalAlignment','center');
        rectangle('Position',boundbox_kk(:,ind),'edgecolor','b','linewidth',1);
        title('Make sure the camera number is right:');
    end;
end;
hold off;