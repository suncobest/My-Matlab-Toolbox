% click_track_mirror: Click the button "Extract grid corners" of the Multimirrors Calibration Toolbox to run this script

if ~exist('strnum_cell','var'),
    read_images_mirror;
end;

fprintf(1,['\n--------------------------------Multimirrors ' ...
           'Initialization--------------------------------\n']);
if ~exist('n_cam','var'),
    n_cam = input('Number of subviews in each image = ([] = 1)');
    if isempty(n_cam),
        n_cam = 1;
    end;
end;

% active_imgviews(i,j) denote the availability of the j-th view of i-th image:1 for yes, 0 for no
active_imgviews = ones(n_cam,1) * active_images;
n_view = n_ima * n_cam;

% hand_list(i) means axis handness of the i-th view
if ~exist('hand_list','var'),
    flag = 1;
    fprintf(1,'\nHandedness vector of camera views do not exist!\n');
else
    hand_list =  hand_list(:)';
    flag = ~isequal(abs(hand_list),ones(1,n_cam));
    if flag,
        fprintf(1,'\nUnexpected value for handedness vector! Please input again!\n');
    end;
end;
while flag,
    fprintf(1,['\nPlease input the axis handedness vector of all ' num2str(n_cam) ' views!\n' ...
                    'Note: 1 stands for right-handed, -1 stands for left-handed.\n']);
    hand_list = input(['handedness vector = ([] = [' num2str(ones(1,n_cam)) '])']);
    if isempty(hand_list),
        hand_list = ones(1,n_cam);
        flag = 0;
    else
        hand_list =  hand_list(:)';
        flag = ~isequal(abs(hand_list),ones(1,n_cam));
        if flag,
            fprintf(1,'\nUnexpected value for handedness vector! Please input again!\n');
        end;
    end;
end;

if ~exist('map','var'),
    map = gray(256);
end;

if ~exist('wintx_default','var') || ~exist('winty_default','var'),
    if ~exist('nx','var'),
        wintx_default = 5;
    else
        wintx_default = max(round(nx/128),round(ny/96));
    end;
    winty_default = wintx_default;
    clear wintx winty;
end;

if ~exist('wintx','var') || ~exist('winty','var'),
    clear win_size; % Clear all the window sizes (to re-initiate)
end;

if ~exist('is_movie','var'),
    fprintf(1,'\nMotion picture (like movie) can use the tracking function:\n');
    is_movie = input('Do you have successive images for corner extraction? ([]=no, other=yes) ','s');
    is_movie = ~isempty(is_movie);
end;

if is_movie,
    same_pattern = 1;
    reset_ask = 0;
    ima_proc = 1:n_ima;
else
    if ~exist('same_pattern','var'),
        same_pattern = input('\nDo you have the same pattern in every image? ([]=yes, other=no) ','s');
        same_pattern = isempty(same_pattern);
    end;
    reset_ask = input('Do you want to reset the image numbers for calibration each time? ([]=no, other=yes) ','s');
    reset_ask = ~isempty(reset_ask);
end;

if same_pattern,
    if exist('dX','var'),
        dX_default = dX;
    end;
    if exist('dY','var'),
        dY_default = dY;
    end;
    if exist('n_sq_x','var'),
        n_sq_x_default = n_sq_x;
    end;
    if exist('n_sq_y','var'),
        n_sq_y_default = n_sq_y;
    end;
end;

if ~exist('dX_default','var') || ~exist('dY_default','var');
    % Setup of JY - 3D calibration rig at Intel (new at Intel) - use units in mm to match Zhang
    dX_default = 30;
    dY_default = 30;

    % Setup of JY - 3D calibration rig at Google - use units in mm to match Zhang
    %     dX_default = 100;
    %     dY_default = 100;
end;

if ~exist('n_sq_x_default','var') || ~exist('n_sq_y_default','var'),
    n_sq_x_default = 10;
    n_sq_y_default = 10;
end;

if reset_ask,
    ima_proc = input('Number(s) of image(s) to process = ([] = all images)');
    if isempty(ima_proc),
        ima_proc = 1:n_ima;
    end;
else
    if ~exist('ima_proc','var'),
        ima_proc = input('Number(s) of image(s) to process = ([] = all images)');
        if isempty(ima_proc),
            ima_proc = 1:n_ima;
        end;
    end;
end;

if exist('win_size', 'var') && ~isequal(size(win_size),[2, n_cam]),
    fprintf(1,'\nUnexpected dimension of window-size list! The window-size list will be reinitiated.\n');
    clear win_size;
end;
if ~exist('win_size', 'var'),
    win_size = zeros(2,n_cam);
    for pp=1:n_cam,
        fprintf(1,'\nWindow size for corner finder (width and height) of view %d:\n', pp);
        wintx = input(['wintx = ([] = ' num2str(wintx_default) ')']);
        if isempty(wintx), wintx = wintx_default; end;
        wintx = round(wintx);
        winty = input(['winty = ([] = ' num2str(winty_default) ')']);
        if isempty(winty), winty = winty_default; end;
        winty = round(winty);
        win_size(:, pp) = [wintx; winty];
        fprintf(1,'Window size of view %d = %dx%d\n', pp,2*wintx+1,2*winty+1);
    end;
end;

if is_movie,
    manual_squares = 0;
    if ~exist('winx_factor','var') || ~exist('winy_factor','var'),
        fprintf(1,['\nWindow size for corner tracker is assumed to be multiple of win_size!\n ' ...
                   'Please input the scaling factor of tracking window:\n']);
        winx_factor = input('winx_factor = ([] = 3)');
        if isempty(winx_factor),
            winx_factor =  3;
        end;
        winy_factor = input('winy_factor = ([] = 3)');
        if isempty(winy_factor),
            winy_factor =  3;
        end;
        winxtrack = wintx * winx_factor;
        winytrack = winty * winy_factor;
        fprintf(1,'Window size for tracking is about %dx%d\n',2*winxtrack+1,2*winytrack+1);
    end;
else
    fprintf(1,['Do you want to use the automatic square counting function\n' ...
        'or do you want to enter the square numbers manually ...\n']);
    manual_squares = input('Enter square numbers manually? ([]=no, other=yes) ','s');
    manual_squares = ~isempty(manual_squares);
end;

if exist('winx_factor','var')~=1,
    winx_factor =  3;
    winy_factor = 3;
end;
fprintf(1,'\nExtraction of the grid corners on the images:\n');

% check if there is extracted corner data
keep_corner = 0;
if  exist('x_cell', 'var') && isequal(size(x_cell),[1, n_view]) && sum(OK_list),
    fprintf(1,['\nExtracted corner data detected!\nDo you want to continue extraction from the ' ...
        'last position or restart from the begining ...\n']);
    keep_corner = input('Continue extraction or not? ([]=continue; other=restart.) ','s');
    keep_corner = isempty(keep_corner);
end;

if ~keep_corner,
    clear x_cell;
end;

if ~exist('x_cell', 'var'),
    x_cell = cell(1,n_view);
    X_cell = x_cell;
    OK_list = false(1,n_view);
    dXY_mat = NaN(2,n_view);
    n_sq_mat = dXY_mat;
    delta = max(nx, ny)/40;
    no_dXY = 1;
end;

for pp = 1:n_cam,
    track_flag = 0;   % indicate tracking or not
    wintx = win_size(1,pp);
    winty = win_size(2,pp);
    winxtrack = wintx * winx_factor;
    winytrack = winty * winy_factor;
    fprintf(1,['Window size of view %d = %dx%d      (Note: To reset the window size, clear the ' ...
        'variable ''win_size'')\n'], pp,2*wintx+1,2*winty+1);

    for kk = ima_proc,
        kth = (kk-1)*n_cam+pp;
        if keep_corner,     % continue extraction or not
            if OK_list(kth)
                % fprintf(1,'\nView %d of image %d is done.\n',pp,kk);
                continue;
            else
                x_cell{kth} = [];
            end;
        end;
        string_num = strnum_cell{kk};    % strnum_cell contain sorted string numbers of images
        ima_name = [calib_name  string_num '.' format_image];

        if exist(ima_name,'file')==2,
            fprintf(1,'\nProcessing with view %d of image %d...\n',pp,kk);
            fprintf(1,'Loading image %s...\n',ima_name);

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

            active_view = 1;
            grid_success = 0;
            while (~grid_success),
                if is_movie && track_flag,
                    manual_click = 0;
                    [Xc,good,outspeed] = cornertracker(Xc0,I0,I,wintx,winty,winxtrack,winytrack,inspeed);
                    if ~all(good),                    % reclick four extreme corners when tracking is bad
                        track_flag = 0;
                    end;
                end;

                if ~track_flag,
                    manual_click = 1;
                    figure(2); clf;
                    image(I);
                    colormap(map);
                    axis image;
                    set(2,'color',[1 1 1]);
                    title(['Check the availability of the view ' num2str(pp) ' of Image ' num2str(kk)]);
                    fprintf(1,['\nIs the view ' num2str(pp) ' of the image ' num2str(kk) ' available ' ...
                                        'for calibration?\nNote: yes=1, no=0. (default=1)\n']);

                    active_view = input(['The availability of the view ' num2str(pp) ' of the image: ([]=1)']);   % check the availability of each view
                    if isempty(active_view)
                        active_view = 1;
                    else
                        active_view = ~~active_view;
                    end;

                    if active_view,
                        bad_clicks = 1;
                        while bad_clicks,
                            figure(2); clf;
                            image(I);
                            colormap(map);
                            axis image;
                            set(2,'color',[1 1 1]);
                            title(['View ' num2str(pp) ' of Image ' num2str(kk) ': Click on the ' ...
                                                'four extreme corners(first corner = origin)...']);
                            disp(['Click on the four extreme corners of the rectangular complete ' ...
                                  'pattern (the first clicked corner is the origin)...']);

                            x= [];y = [];
                            figure(2); hold on;
                            for count = 1:4,
                                [xi,yi] = ginput(1);
                                [xxi] = cornerfinder([xi;yi],I,wintx,winty);
                                xi = xxi(1);
                                yi = xxi(2);
                                figure(2);
                                plot(xi,yi,'+','color',[ 1.000 0.314 0.510 ],'linewidth',2);
                                plot(xi+[wintx+.5, -(wintx+.5), -(wintx+.5), wintx+.5, wintx+.5], ...
                                      yi+[winty+.5, winty+.5, -(winty+.5), -(winty+.5), winty+.5], ...
                                       '-','color',[ 1.000 0.314 0.510 ],'linewidth',2);
                                x = [x;xi];
                                y = [y;yi];
                                plot(x,y,'-','color',[ 1.000 0.314 0.510 ],'linewidth',2);
                                drawnow;
                            end;
                            plot([x;x(1)],[y;y(1)],'-','color',[ 1.000 0.314 0.510 ],'linewidth',2);
                            drawnow;
                            hold off;
                            [Xc,good,bad,type] = cornerfinder([x';y'],I,wintx,winty); % the four corners
                            bad_clicks = (sum(bad)>0);
                        end;
                    else
                        grid_success = 1;            % get out of the loop
                    end;
                end;

                if active_view,
                    x = Xc(1,:)';
                    y = Xc(2,:)';

                    % Sort the corners:
                    x_mean = mean(x);
                    y_mean = mean(y);
                    x_v = x - x_mean;
                    y_v = y - y_mean;
                    % atan(-y,x)=-atan(y,x); handedness: right=1,left=-1
                    theta = -hand_list(pp) * atan2(y_v,x_v);
                    % ��theta = atan2(y_v,x_v)��Ӧ�ķ���Ϊ(a1,a2,a3,a4), �õ��ǶȲ�theta-theta(1)Ϊ
                    % (a1-a1,a2-a1,a3-a1,a4-a1)=angle,��theta = atan2(-y_v,x_v)�ķ��Ƕ�Ϊ
                    % (-a1,-a2,-a3,-a4)���ǶȲ�theta-theta(1)Ϊ(a1-a1,a1-a2,a1-a3,a1-a4)=-angle. ��
                    % ͼ������ϵ�У����ǵ�������(��x����ת��ʸ��)Ϊ˳ʱ��angle�����߳��������ͼ��
                    % ����������ϵ��xyƽ����ͼ���ϵ�ͶӰ��ʱ��Ϊ��(-angle)�����߳��⡣
                    % mod(x,y)=x-n.*y ���� n=floor(x./y)�����ڽǶȲ�a����-2*pi<=a<=2*pi����a>=0,��
                    % mod(a,2*pi)=a����a<0,��mod(a,2*pi)=2*pi+a������0<=mod(a,2*pi)<2*pi����a������
                    % [0,2*pi)����-angle���������Ǵ�ԭ����ʱ������,angle��֮�������-angle�����
                    % ����ʱ��Ϊ���������ϵ(handedness=1)����angle����õ�˳ʱ��Ϊ���������ϵ(hand=-1)��

                    [junk,ind] = sort(mod(theta-theta(1),2*pi));
                    % ind: counterclockwise ascending order; first point is the origin (o,x,xy,y)
                    ind = ind([4 3 2 1]);
                    % ind: descending order; last point is the origin (y,xy,x,o)
                    x = x(ind);
                    y = y(ind);
                    % (x1,y1) is the extreme point in y direction, (x3,y3) is the extreme point in x direction
                    % (x4,y4) is the origin, (x2,y2) is the extreme point in the diagonal direction
                    x1= x(1); x2 = x(2); x3 = x(3); x4 = x(4);
                    y1= y(1); y2 = y(2); y3 = y(3); y4 = y(4);

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

                    figure(2); image(I);
                    axis image;
                    colormap(map);
                    hold on;
                    plot([x;x(1)],[y;y(1)],'g-');
                    plot(x,y,'og');
                    hx=text(x6 + delta * vX(1) ,y6 + delta*vX(2),'X');
                    set(hx,'color','g','Fontsize',14,'HorizontalAlignment','center');
                    hy=text(x7 + delta*vY(1), y7 + delta*vY(2),'Y');
                    set(hy,'color','g','Fontsize',14,'HorizontalAlignment','center');
                    hO=text(x4 + delta * vO(1) ,y4 + delta*vO(2),'O','color','g','Fontsize',14, ...
                            'HorizontalAlignment','center');
                    for ii = 1:4,
                        text(x(ii),y(ii),num2str(ii));
                    end;
                    hold off;

                    if manual_squares,
                        n_sq_x = input(['Number of squares along the X direction = ([]=' num2str(n_sq_x_default) ')']);
                        if isempty(n_sq_x), n_sq_x = n_sq_x_default; end;
                        n_sq_y = input(['Number of squares along the Y direction = ([]=' num2str(n_sq_y_default) ')']);
                        if isempty(n_sq_y), n_sq_y = n_sq_y_default; end;
                        grid_success = 1;
                    else
                        % Try to automatically count the number of squares in the grid
                        n_sq_x1 = count_squares(I,x1,y1,x2,y2,wintx);
                        n_sq_x2 = count_squares(I,x3,y3,x4,y4,wintx);
                        n_sq_y1 = count_squares(I,x2,y2,x3,y3,wintx);
                        n_sq_y2 = count_squares(I,x4,y4,x1,y1,wintx);

                        % If could not count the number of squares, enter manually
                        if n_sq_x1~=n_sq_x2 || n_sq_y1~=n_sq_y2 || min([n_sq_x1 n_sq_x2 n_sq_y1 n_sq_y2])<0,
                            disp('Could not count the number of squares in the grid.');
                            flag = input('Do you want to re-click the four corners of this view? ([]=1=yes) >>');
                            if isempty(flag) || flag,
                                grid_success = 0;
                            else
                                disp('Enter the number of squares manually.');
                                n_sq_x = input(['Number of squares along the X direction = ([]=' num2str(n_sq_x_default) ')']);
                                if isempty(n_sq_x), n_sq_x = n_sq_x_default; end;
                                n_sq_y = input(['Number of squares along the Y direction = ([]=' num2str(n_sq_y_default) ')']);
                                if isempty(n_sq_y), n_sq_y = n_sq_y_default; end;
                                grid_success = 1;
                            end;
                        else
                            n_sq_x = n_sq_x1;
                            n_sq_y = n_sq_y1;
                            grid_success = 1;
                        end;
                    end;

                    if is_movie,
                        if grid_success,
                            I0 = I;
                            Xc0 = Xc;
                            if track_flag,
                                inspeed = outspeed;                % good tracking: update speed
                            else
                                inspeed = zeros(2,4);              % manual_click: reset speed to 0
                            end;
                            track_flag = 1;
                            fprintf(1,'Switch tracking function on!\n');
                        else
                            track_flag = 0;
                            fprintf(1,'\nBad tracking! Switch tracking function off!\n\n');
                        end;
                    end;

                    if ~grid_success
                        fprintf(1,'Invalid grid. Try again.\n');
                    end;
                end;
            end;

            if active_view,
                n_sq_x_default = n_sq_x;
                n_sq_y_default = n_sq_y;
                if ~same_pattern || no_dXY,
                    % Enter the size of each square
                    dX = input(['Size dX of each square along the X direction = ([]=' num2str(dX_default) 'mm)']);
                    dY = input(['Size dY of each square along the Y direction = ([]=' num2str(dY_default) 'mm)']);
                    if isempty(dX), dX = dX_default; else dX_default = dX; end;
                    if isempty(dY), dY = dY_default; else dY_default = dY; end;
                    no_dXY = 0;
                else
                    fprintf(1,['Size of each square along the X direction: dX=' num2str(dX) 'mm\n']);
                    fprintf(1,['Size of each square along the Y direction: dY=' num2str(dY) 'mm\n'...
                               '(Note: To reset the size of the squares, run the script ''no_dXY=1;'')\n']);
                end;

                % Compute the inside points through computation of the planar homography (collineation)
                a00 = [x(1);y(1);1];  % y
                a10 = [x(2);y(2);1];  % xy
                a11 = [x(3);y(3);1];  % x
                a01 = [x(4);y(4);1];  % o

                % Compute the planar collineation:
                Homo = compute_homography_lm([a00 a10 a11 a01],[0 1 1 0;0 0 1 1;1 1 1 1]);

                % Build the grid using the planar collineation:
                x_l = ((0:n_sq_x)'*ones(1,n_sq_y+1))/n_sq_x;
                y_l = (ones(n_sq_x+1,1)*(0:n_sq_y))/n_sq_y;
                % pts array: row number increase along x direction, column number increase along y derection,
                % array scan the x direction: 1st points is y (on the top right corner), origin is on the top left corner.
                pts = [x_l(:) y_l(:) ones((n_sq_x+1)*(n_sq_y+1),1)]';
                XX = Homo*pts;
                XX = XX(1:2,:) ./ (ones(2,1)*XX(3,:));


                %%%%%%%%%%%%%%%% ADDITIONAL STUFF IN THE CASE OF HIGHLY DISTORTED IMAGES %%%%%%%%%%%%%
                figure(2);
                hold on;
                plot(XX(1,:),XX(2,:),'r+');
                title('The red crosses should be close to the image corners');
                hold off;

                % Complete size of the rectangle
                W = n_sq_x*dX;
                L = n_sq_y*dY;
                grid_pts = [x(1) x(2) x(4) x(3);y(1) y(2) y(4) y(3)];           % [y xy o x]
                Xgrid = [0 W 0 W;L L 0 0];

                if manual_click,     %% guess distortion
                    disp('If the guessed grid corners (red crosses on the image) are not close to the actual corners,');
                    disp(['it is necessary to enter an initial guess for the radial distortion ' ...
                          'factor kc (useful for subpixel detection)']);
                    quest_distort = input('Need of an initial guess for distortion? ([]=no, other=yes) ','s');
                    quest_distort = ~isempty(quest_distort);
                    if quest_distort,
                        % Estimation of focal length:
                        c_g = [nx;ny]/2 + .5;
                        f_g = focus_extrinsic_estimate(grid_pts,Xgrid,c_g);
                        f_g = mean(f_g);

                        satis_distort = 0;
                        disp(['Estimated focal: ' num2str(f_g) ' pixels']);

                        while ~satis_distort,
                            k_g = input('Guess for distortion factor kc ([]=0): ');
                            if isempty(k_g), k_g = 0; end;
                            xy_corners_undist = comp_distortion2([x' - c_g(1);y'-c_g(2)]/f_g,k_g);
                            xu = xy_corners_undist(1,:)';
                            yu = xy_corners_undist(2,:)';
                            [XXu] = projectedGrid ( [xu(1);yu(1)], [xu(2);yu(2)],[xu(3);yu(3)], ...
                                                    [xu(4);yu(4)],n_sq_x+1,n_sq_y+1); % The full grid
                            XX = (ones(2,1)*(1 + k_g * sum(XXu.^2))) .* XXu;
                            XX(1,:) = f_g*XX(1,:)+c_g(1);
                            XX(2,:) = f_g*XX(2,:)+c_g(2);

                            figure(2);
                            image(I);
                            colormap(map);
                            zoom on;
                            hold on;
                            plot(XX(1,:),XX(2,:),'r+');
                            title('The red crosses should be on the grid corners...');
                            hold off;
                            satis_distort = input('Satisfied with distortion? ([]=no, other=yes) ','s');
                            satis_distort = ~isempty(satis_distort);
                        end;
                    end;

                else          % do not guess distortion, leave some parameters unchanged (quest_distort,f_g,c_g,k_g)

                    if quest_distort,
                        f_g = focus_extrinsic_estimate(grid_pts,Xgrid,c_g,f_g,k_g,0);
                        f_g = mean(f_g);
                        xy_corners_undist = comp_distortion2([x' - c_g(1);y'-c_g(2)]/f_g,k_g);              % highly accurate
                        xu = xy_corners_undist(1,:)';
                        yu = xy_corners_undist(2,:)';

                        [XXu] = projectedGrid ( [xu(1);yu(1)], [xu(2);yu(2)],[xu(3);yu(3)], ...
                                                [xu(4);yu(4)],n_sq_x+1,n_sq_y+1); % The full grid

                        XX = (ones(2,1)*(1 + k_g * sum(XXu.^2))) .* XXu;
                        XX(1,:) = f_g*XX(1,:)+c_g(1);
                        XX(2,:) = f_g*XX(2,:)+c_g(2);

                        figure(2);
                        image(I);
                        colormap(map);
                        zoom on;
                        hold on;
                        plot(XX(1,:),XX(2,:),'r+');
                        title('The red crosses should be on the grid corners...');
                        hold off;
                    end;
                end;
                %%%%%%%%%%%%% END ADDITIONAL STUFF IN THE CASE OF HIGHLY DISTORTED IMAGES %%%%%%%%%%%%%

                Np = (n_sq_x+1)*(n_sq_y+1);
                disp('Corner extraction...');
                %%% Finds the exact corners at every points!
                % subtract 1 to bring the origin to (0,0) instead of (1,1) in matlab (not necessary in C)
                grid_pts = cornerfinder(XX,I,wintx,winty);
                grid_pts = grid_pts - 1;

                ind_corners = [1, n_sq_x+1, (n_sq_x+1)*n_sq_y+1, (n_sq_x+1)*(n_sq_y+1)];
                % index of the 4 corners(y,xy,o,x)
                ind_orig = (n_sq_x+1)*n_sq_y + 1;
                xorig = grid_pts(1,ind_orig);
                yorig = grid_pts(2,ind_orig);
                dxpos = mean([grid_pts(:,ind_orig) grid_pts(:,ind_orig+1)],2);
                dypos = mean([grid_pts(:,ind_orig) grid_pts(:,ind_orig-n_sq_x-1)],2);

                x_box_kk = [grid_pts(1,:)-(wintx+.5); grid_pts(1,:)+(wintx+.5); grid_pts(1,:)+(wintx+.5); ...
                    grid_pts(1,:)-(wintx+.5); grid_pts(1,:)-(wintx+.5)];
                y_box_kk = [grid_pts(2,:)-(winty+.5); grid_pts(2,:)-(winty+.5); grid_pts(2,:)+(winty+.5); ...
                    grid_pts(2,:)+(winty+.5); grid_pts(2,:)-(winty+.5)];

                figure(3);
                image(I); axis image; colormap(map); hold on;
                plot(grid_pts(1,:)+1,grid_pts(2,:)+1,'r+');
                plot(x_box_kk+1,y_box_kk+1,'-b');
                plot(grid_pts(1,ind_corners)+1,grid_pts(2,ind_corners)+1,'mo');
                plot(xorig+1,yorig+1,'*m');
                h = text(xorig+delta*vO(1),yorig+delta*vO(2),'O');
                set(h,'Color','g','FontSize',14,'HorizontalAlignment','center');
                h2 = text(dxpos(1)+delta*vX(1),dxpos(2)+delta*vX(2),'dX');
                set(h2,'Color','g','FontSize',14,'HorizontalAlignment','center');
                h3 = text(dypos(1)+delta*vY(1),dypos(2)+delta*vY(2),'dY');
                set(h3,'Color','g','FontSize',14,'HorizontalAlignment','center');
                title(['Extracted corners in view ' num2str(pp) ' of image ' num2str(kk)]);
                zoom on;
                drawnow;
                hold off;

                %  the 3D world coordinates
                Xi = reshape(((0:n_sq_x)*dX)'*ones(1,n_sq_y+1),1,Np);
                Yi = reshape(ones(n_sq_x+1,1)*(n_sq_y:-1:0)*dY,1,Np);

                % All the point coordinates (on the image, and in 3D) - for global optimization:
                x = grid_pts;
                X = [Xi;Yi];

                % Saves all the data into variables:
                x_cell{kth}= x;
                X_cell{kth}= X;

                dXY_mat(1, kth)  = dX;
                dXY_mat(2, kth)  = dY;

                n_sq_mat(1,kth) = n_sq_x;
                n_sq_mat(2,kth)  = n_sq_y;
            end;

        else  % if the kk-th image do not exist, then variable keep untouched: [] for cell, 'NaN' for mat
            active_view = 0;
            if is_movie,
                track_flag = 0;
                fprintf(1,'\nImage don''t exist! Switch tracking function off!\n\n');
            end;
        end;
        active_imgviews(pp,kk) = active_view;
        OK_list(kth) = 1;
    end;
end;

check_extracted_mirror;

string_save = ['save multimirror_calib_data active_images ind_active wintx winty n_ima format_image ' ...
               'calib_name nx ny map dX_default dY_default dX dY X_cell x_cell n_sq_mat dXY_mat ' ...
               'OK_list win_size n_cam active_imgviews hand_list strnum_cell'];
eval(string_save);
disp('done.');
fprintf(1,'NOTE: Run script ''swap_corner_xy'' if your (transparent) checkerboard turn upside down in one view!\n');
