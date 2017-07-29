% recheck_corners_mirror
% This script will draw corners on every image, please check the corners and world axes.
% Make sure that world axes in each view of each image are consistent with each other.

if ~exist('strnum_cell','var'),
    if exist('multimirror_calib_data.mat','file'),
        load('multimirror_calib_data.mat');
    else
        fprintf(1,['\nThere is no data available to recheck corners! \nPlease run "multimirrors_gui" '...
            'and click the first two buttons!\n"Read images" , then "Extract grid corners"...\n']);
        return;
    end;
end;

check_extracted_mirror;

for kk = ind_active,
    string_num = strnum_cell{kk};                           % strnum_cellÎª°´Ë³ÐòÅÅÁÐµÄÍ¼Æ¬ÐòºÅ×Ö·ûcell
    ima_name = [calib_name  string_num '.' format_image];

    if ~exist(ima_name,'file'),
        fprintf(1,'\nWarning: Image %s not found in the directory!\n',ima_name);
        return;
    end;
end;

if  length(ind_active)>30,
    flag = input('Do you want to check images one by one or play the entire sequence? ([]=manual, other=auto)' ,'s');
    flag = ~isempty(flag);
else
    flag = 0;
end;
delta = max(nx, ny)/40;        % textÆ«ÒÆÁ¿
% set waiting time to refresh images
waitDt = 0.3;

for kk = ind_active,
    string_num = strnum_cell{kk};                           % strnum_cellÎª°´Ë³ÐòÅÅÁÐµÄÍ¼Æ¬ÐòºÅ×Ö·ûcell
    ima_name = [calib_name  string_num '.' format_image];

    if format_image(1) == 'p',
        if format_image(2) == 'p',
            I = double(loadppm(ima_name));
        elseif format_image(2) == 'g',
            I = double(loadpgm(ima_name));
        else
            I = double(imread(ima_name));
        end;
    else
        if format_image(1) == 'r',
            I = readras(ima_name);
        else
            I = double(imread(ima_name));
        end;
    end;
    if size(I,3)==3,
        I = 0.299 * I(:,:,1) + 0.5870 * I(:,:,2) + 0.114 * I(:,:,3);
    end;

    figure(2);
    image(I); axis image; colormap(gray(256)); hold on;
    title(['Check the extracted corners on image: ' num2str(kk)]);
    active_view = active_imgviews(:,kk);
    for pp = 1:n_cam,
        if active_view(pp),
            kth = (kk-1)*n_cam+pp;
            x_kk = x_cell{kth};
            n_sq_x = n_sq_mat(1,kth);
            n_sq_y = n_sq_mat(2,kth);

            ind_corners = [1, n_sq_x+1, (n_sq_x+1)*(n_sq_y+1), (n_sq_x+1)*n_sq_y+1];   % index of the 4 corners(y,xy,x,o)
            x1 = x_kk(1,1);                     % y
            y1 = x_kk(2,1);
            x2 = x_kk(1,ind_corners(2));        % xy
            y2 = x_kk(2,ind_corners(2));
            x3 = x_kk(1,ind_corners(3));        % x
            y3 = x_kk(2,ind_corners(3));
            x4 = x_kk(1,ind_corners(4));        % o
            y4 = x_kk(2,ind_corners(4));

            % Find center: (x5,y5)
            p_center = cross(cross([x1;y1;1],[x3;y3;1]),cross([x2;y2;1],[x4;y4;1]));
            x5 = p_center(1)/p_center(3);
            y5 = p_center(2)/p_center(3);
            % center on the X axis:
            x6 = (x3 + x4)/2;
            y6 = (y3 + y4)/2;
            % center on the Y axis:
            x7 = (x1 + x4)/2;
            y7 = (y1 + y4)/2;
            % Direction of displacement for the X axis:
            vX = [x6-x5;y6-y5];
            vX = vX / norm(vX);
            % Direction of displacement for the Y axis:
            vY = [x7-x5;y7-y5];
            vY = vY / norm(vY);
            % Direction of diagonal:
            vO = [x4 - x5; y4 - y5];
            vO = vO / norm(vO);

            % plot corners on each view
            fprintf(1,'\nPloting extracted corners on view %d of image %d...\n',pp,kk);
            figure(2);
            plot(x_kk(1,:)+1,x_kk(2,:)+1,'r+');
            plot(x_kk(1,[ind_corners,1])+1, x_kk(2,[ind_corners,1])+1,'g-');
            plot(x_kk(1,ind_corners)+1, x_kk(2,ind_corners)+1,'mo');
            plot(x4+1,y4+1,'m*');
            h = text(x4+delta*vO(1),y4+delta*vO(2),'O');
            set(h,'Color','g','FontSize',14,'HorizontalAlignment','center');
            h2 = text(x6+delta*vX(1),y6+delta*vX(2),'X');
            set(h2,'Color','g','FontSize',14,'HorizontalAlignment','center');
            h3 = text(x7+delta*vY(1),y7+delta*vY(2),'Y');
            set(h3,'Color','g','FontSize',14,'HorizontalAlignment','center');
        end;
    end;
    % set figure 2
    set(2,'color',[1 1 1],'Name',num2str(kk),'NumberTitle','off');
    zoom on;
    drawnow;
    hold off;
    if flag,
        pause(waitDt);
    else
        fprintf(1,'\nCheck the extracted corners on image %d, press any key to continue...\n', kk);
        pause;
    end;
end;

fprintf(1, ['\ndone.\n\nNote: if the corner on view ''pp'' of image ''kk'' is wrongly extracted,\n' ...
    'delete the extracted data of the image by running the script ''Ok_list((kk-1)*n_cam+pp) = 0;''\n'...
    'and extract corner again!\n']);
