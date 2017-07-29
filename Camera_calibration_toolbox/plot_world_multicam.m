%%%%%%%%%%% SHOW EXTRINSIC RESULTS  IN WORLD FRAME%%%%%%%%%%%%%

% Qcc (Omcc),Tcc,handcc是所有摄像机相对于第一个摄像机的运动变换，所以第一列
% 为cam1到自己的变换; 在第一个摄像机坐标系中画出所有摄像机坐标系，以及标定板

delta = 0.2;        % shift ot text
if ~exist('n_ima','var')||~exist('fc_mat','var')||~exist('Tw_mat','var'),
    fprintf(1,'No  calibration data available.\n');
    return;
end;

ind_active_views = find(active_imgviews(:)');
% unit: normT为世界坐标系原点到摄像机坐标系原点最大距离的1/10
normT = sqrt(max(sum(Tw_mat(:,ind_active_views).^2,1),[],2))/10;
no_grid = 0;   % rig3D = 1;
if ~exist('n_sq_mat','var'),
    no_grid = 1;
    dX = normT;
    dY = dX;
end;

if ~exist('show_camera','var'),
    show_camera = 1;
end;

if ~exist('refine_multicam','var'),
    if n_cam>1,
        refine_multicam = 1;
    else
        refine_multicam = 0;
    end;
end;

if ~exist('Qw_mat','var') && exist('Omw_mat','var'),
    Qw_mat = trans_quat_axis(Omw_mat);
    if refine_multicam,
        Qcw = trans_quat_axis(Omcw);
        Qcc = trans_quat_axis(Omcc);
    end;
end;

if show_camera,
    if ~exist('hand1','var')||~exist('handcc','var'),
        hand1 = hand_list(1);
        handcc = hand1*hand_list;
    end;
    BASE = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]*normT;  % 原点 x基矢量 原点 y基矢量 原点 z基矢量
    DeltaXYZ = [BASE(:,[2,4,6])*(1+delta), -[1;1;1]*normT*delta];     % text位置：x, y, z, o

    % 求出图片像素坐标的四个顶点(0,0), (nx-1,0), (nx-1,ny-1), (0,ny-1)的摄像机归一化坐标
    % 若xp和xd都是齐次坐标(3*N, 第三行为1), 则有 xp= KK*xd; 所以xd=KK\xp.
    % camera axes base vector initialization (wrt camera 1)
    BASE_cell=cell(1,n_cam);
    DELTA_cell=BASE_cell;
    IP_cell=BASE_cell;
    active_view = find(any(active_imgviews,2)');
    for pp = active_view,
        fc = fc_mat(:,pp);
        cc = cc_mat(:,pp);
        alpha_c = alpha_vec(pp);
        nx = imsize(1,pp);
        ny = imsize(2,pp);
        KK = [fc(1) fc(1)*alpha_c cc(1);0 fc(2) cc(2); 0 0 1];
        IP = KK\[0 nx-1 nx-1 0 0 ; 0 0 ny-1 ny-1 0;1 1 1 1 1]*normT;
        IP = reshape([IP;BASE(:,1)*ones(1,5);IP],3,15);
        % Change of reference: wrt camera1
        Rckk = trans_quat_mat(Qcc(:,pp));
        if handcc(pp)~=1,
            Rckk(:,3)=-Rckk(:,3);
        end;
        Twkk = Tcc(:,pp);
        BASE_cell{pp} = Rckk'*(BASE-Twkk(:,ones(1,6)));       % X1=Rk1'*(Xk-Tk1)
        DELTA_cell{pp}=Rckk'*(DeltaXYZ-Twkk(:,ones(1,4)));
        IP_cell{pp} = Rckk'*(IP -Twkk(:,ones(1,15)));
    end;
end;

figure(4);
clf;
hold on;
% Color code for each view:
palette = 'brgkcm';
if exist('X_cell','var'),
     ind_active_views = find(active_imgviews(1,:));
    for kk = ind_active_views,
        kth =   (kk-1)*n_cam+1;
        Y_kk = X_cell{kth};
        [m,n] = size(Y_kk);
        if m==2,
            Y_kk = [Y_kk; zeros(1,n)];
        end;
        if ~no_grid,
            n_sq_x = n_sq_mat(1,kth);
            n_sq_y = n_sq_mat(2,kth);
            dX = dXY_mat(1,kth);
            dY = dXY_mat(2,kth);
            if (n ~= (n_sq_x+1)*(n_sq_y+1)),
                no_grid = 1;
            end;
        end;
        figure(4);
        if ~no_grid,
            Yx = reshape(Y_kk(1,:),n_sq_x+1,n_sq_y+1);
            Yy = reshape(Y_kk(2,:),n_sq_x+1,n_sq_y+1);
            Yz = reshape(Y_kk(3,:),n_sq_x+1,n_sq_y+1);
            mesh(Yx,Yy,Yz,'edgecolor','k','facecolor','none','linewidth',2);
        else
            plot3(Y_kk(1,:),Y_kk(2,:),Y_kk(3,:),'k.');
        end;
        if show_camera,
            if refine_multicam,
                Qwkk = Qcw(:,kk);
                Twkk = Tcw(:,kk);
            else
                Qwkk = Qw_mat(:,kth);
                Twkk = Tw_mat(:,kth);
            end;
            % Change of reference: wrt world frame
            Rckk = trans_quat_mat(Qwkk);
            if hand1~=1,
                Rckk(:,3)=-Rckk(:,3);
            end;
            % draw camers in world frame
            % 按照Matlab坐标系特点，将摄像机坐标系基矢量变为(e1, e3, -e2),相当于绕e1旋转-pi/2，使深度方向水平
            active_view = find(active_imgviews(:,kk)');
            for pp = active_view,
                kth =   (kk-1)*n_cam+pp;
                rgbi = palette(rem(pp-1,6)+1);
                BASE_kk = Rckk'*(BASE_cell{pp}-Twkk(:,ones(1,6)));       % Xw=Rc1'*(X1-T1)
                IP_kk = Rckk'*(IP_cell{pp} -Twkk(:,ones(1,15)));
                DELTA_kk = Rckk'*(DELTA_cell{pp}-Twkk(:,ones(1,4)));
                plot3(BASE_kk(1,:),BASE_kk(2,:),BASE_kk(3,:),'-','color',rgbi,'linewidth',1);
                plot3(IP_kk(1,:),IP_kk(2,:),IP_kk(3,:),'-','color',rgbi,'linewidth',1);
                % show face
                patch(struct('vertices',IP_kk','faces',[1 2 4;4 2 7;7 2 10;10 2 1]),'facecolor',[52 217 160]/255,'EdgeColor', 'r');
                patch(struct('vertices',IP_kk','faces',[1 4 7 10]),'facecolor',[247 239 7]/255,'EdgeColor', 'none');
                % show text
                text(DELTA_kk(1,1),DELTA_kk(2,1),DELTA_kk(3,1),'X','color',rgbi,'HorizontalAlignment','center');
                text(DELTA_kk(1,2),DELTA_kk(2,2),DELTA_kk(3,2),'Y','color',rgbi,'HorizontalAlignment','center');
                text(DELTA_kk(1,3),DELTA_kk(2,3),DELTA_kk(3,3),'Z','color',rgbi,'HorizontalAlignment','center');
                text(DELTA_kk(1,4),DELTA_kk(2,4),DELTA_kk(3,4),num2str(kth),'color',rgbi,'HorizontalAlignment','center');
            end;
        end;
    end;
end;

az = 50;
el = 20;
figure(4);
rotate3d on;
grid on;
title('Extrinsic results in world frame');
set(4,'color',[1 1 1],'Name','3D','NumberTitle','off');
axis equal vis3d tight;
xlabel('X_{world}');
ylabel('Y_{world}');
zlabel('Z_{world}');
view(az,el);
hold off;

if exist('h_sw1','var') && ishandle(h_sw1),
    delete(h_sw1);
end;
if exist('h_sw2','var') && ishandle(h_sw2),
    delete(h_sw2);
end;
t1 = 0.3;       % position: [left bottom width height]
t2 = 0.04;

h_sw1 = uicontrol('Parent',4,'Units','normalized', 'Callback','plot_camera_multicam;', 'Position',[1-t1,0,t1,t2], ...
                           'String','Switch to camera-centered view','fontsize',8,'fontname','clean','Tag','Pushbutton1');
if show_camera,
    h_sw2 = uicontrol('Parent',4,'Units','normalized', 'Callback','show_camera=0;plot_world_multicam;', 'Position', ...
                   [1-t1,t2,t1,t2],'String','Remove camera reference frames','fontsize',8,'fontname','clean','Tag','Pushbutton1');
else
    h_sw2 = uicontrol('Parent',4,'Units','normalized', 'Callback','show_camera=1;plot_world_multicam;', 'Position', ...
                   [1-t1,t2,t1,t2],'String','Add camera reference frames','fontsize',8,'fontname','clean','Tag','Pushbutton1');
end;
