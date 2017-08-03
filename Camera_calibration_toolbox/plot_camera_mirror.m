%%%%%%%%%%% SHOW EXTRINSIC RESULTS  IN VIEW 1 FRAME%%%%%%%%%%%%%

% Qcc (Omcc),Tcc,handcc是所有摄像机相对于第一个摄像机的运动变换，所以第一列
% 为cam1到自己的变换; 在第一个摄像机坐标系中画出所有摄像机坐标系，以及标定板

delta = 0.2;        % shift ot text
if ~exist('n_ima','var')||~exist('fc','var')||~exist('Tw_mat','var'),
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

figure(4);
clf;
hold on;
if show_camera,
    if ~exist('hand1','var')||~exist('handcc','var'),
        hand1 = hand_list(1);
        handcc = hand1*hand_list;
    end;
    if ~exist('KK','var'),
        KK = [fc(1) fc(1)*alpha_c cc(1);0 fc(2) cc(2); 0 0 1];
    end;
    % 求出图片像素坐标的四个顶点(0,0), (nx-1,0), (nx-1,ny-1), (0,ny-1)的摄像机归一化坐标
    % 若xp和xd都是齐次坐标(3*N, 第三行为1), 则有 xp= KK*xd; 所以xd=KK\xp.
    BASE = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]*normT;  % 原点 x基矢量 原点 y基矢量 原点 z基矢量
    DeltaXYZ = [BASE(:,[2,4,6])*(1+delta), -[1;1;1]*normT*delta];     % text位置: x, y, z, o
    IP = KK\[0 nx-1 nx-1 0 0 ; 0 0 ny-1 ny-1 0;1 1 1 1 1]*normT;
    IP = reshape([IP;BASE(:,1)*ones(1,5);IP],3,15);                 % 视线椎

    % camera axes base vector initialization (wrt camera 1)
    % draw camers
    % 按照Matlab坐标系特点，将摄像机坐标系基矢量变为(e1, e3, -e2),相当于绕e1旋转-pi/2，使深度方向水平
    active_view = find(any(active_imgviews,2)');
    for pp = active_view,      % draw camers
        % Change of reference: wrt camera1
        Rckk = trans_quat_mat(Qcc(:,pp));
        if handcc(pp)~=1,
            Rckk(:,3)=-Rckk(:,3);
        end;
        Twkk = Tcc(:,pp);
        BASE_kk = Rckk'*(BASE-Twkk(:,ones(1,6)));       % X1=Rk1'*(Xk-Tk1)
        DELTA_kk = Rckk'*(DeltaXYZ-Twkk(:,ones(1,4)));
        IP_kk = Rckk'*(IP -Twkk(:,ones(1,15)));
        figure(4);
        plot3(BASE_kk(1,:),BASE_kk(3,:),-BASE_kk(2,:),'b-','linewidth',2);
        plot3(IP_kk(1,:),IP_kk(3,:),-IP_kk(2,:),'r-','linewidth',2);
        text(DELTA_kk(1,1),DELTA_kk(3,1),-DELTA_kk(2,1),'X','HorizontalAlignment','center','FontWeight','bold');
        text(DELTA_kk(1,2),DELTA_kk(3,2),-DELTA_kk(2,2),'Y','HorizontalAlignment','center','FontWeight','bold');
        text(DELTA_kk(1,3),DELTA_kk(3,3),-DELTA_kk(2,3),'Z','HorizontalAlignment','center','FontWeight','bold');
        text(DELTA_kk(1,4),DELTA_kk(3,4),-DELTA_kk(2,4),['Cam_' num2str(pp)],'HorizontalAlignment','center','FontWeight','bold');
    end;
end;

% draw mesh or 3D pts
% Color code for each image:
palette = 'brgkcm';

if exist('X_cell','var'),
    ind_active_views = find(active_imgviews(1,:));
    for kk = ind_active_views,
        kth =   (kk-1)*n_cam+1;
        if refine_multicam,
            Qwkk = Qcw(:,kk);
            Twkk = Tcw(:,kk);
        else
            Qwkk = Qw_mat(:,kth);
            Twkk = Tw_mat(:,kth);
        end;
        X_kk = X_cell{kth};
        N_kk = size(X_kk,2);
        if ~no_grid,
            n_sq_x = n_sq_mat(1,kth);
            n_sq_y = n_sq_mat(2,kth);
            dX = dXY_mat(1,kth);
            dY = dXY_mat(2,kth);
            if (N_kk ~= (n_sq_x+1)*(n_sq_y+1)),
                no_grid = 1;
            end;
        end;
        Y_kk = rigid_trans(X_kk,Qwkk,Twkk,hand1);
        uu = [-dX;-dY;0]/2;    % denote orgin of each mesh
        uu = rigid_trans(uu,Qwkk,Twkk,hand1);
        rgbi = palette(rem(kk-1,6)+1);
        figure(4);
        if ~no_grid,
            Yx = reshape(Y_kk(1,:),n_sq_x+1,n_sq_y+1);
            Yy = reshape(Y_kk(2,:),n_sq_x+1,n_sq_y+1);
            Yz = reshape(Y_kk(3,:),n_sq_x+1,n_sq_y+1);
            mesh(Yx,Yz,-Yy,'edgecolor',rgbi,'facecolor','none','linewidth',1);
            text(uu(1),uu(3),-uu(2),num2str(kk),'fontsize',14,'color',rgbi,'HorizontalAlignment','center');
        else
            plot3(Y_kk(1,:),Y_kk(3,:),-Y_kk(2,:),[rgbi '.']);
            text(uu(1),uu(3),-uu(2),num2str(kk),'fontsize',14,'color',rgbi,'HorizontalAlignment','center');
        end;
    end;
end;

az = 50;
el = 20;
figure(4);
rotate3d on;
grid on;
title('Extrinsic results in view 1');
set(4,'color',[1 1 1],'Name','3D','NumberTitle','off');
axis equal vis3d tight;
xlabel('X_{cam1}');
ylabel('Z_{cam1}');
zlabel('-Y_{cam1}');
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

h_sw1 = uicontrol('Parent',4,'Units','normalized', 'Callback','plot_world_mirror;', 'Position',[1-t1,0,t1,t2], ...
                           'String','Switch to world-centered view','fontsize',8,'fontname','clean','Tag','Pushbutton1');
if show_camera,
    h_sw2 = uicontrol('Parent',4,'Units','normalized', 'Callback','show_camera=0;plot_camera_mirror;','Position',...
                   [1-t1,t2,t1,t2],'String','Remove camera reference frames','fontsize',8,'fontname','clean','Tag','Pushbutton1');
else
    h_sw2 = uicontrol('Parent',4,'Units','normalized', 'Callback','show_camera=1;plot_camera_mirror;','Position',...
                   [1-t1,t2,t1,t2],'String','Add camera reference frames','fontsize',8,'fontname','clean','Tag','Pushbutton1');
end;
