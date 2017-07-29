% Omcw,Tcw,hand_list是所有摄像机相对于参考系的运动变换

delta = 0.2;        % shift ot text
if ~exist('n_ima','var') || ~exist('fc_mat','var') || ~exist('Omcw','var'),
    fprintf(1,'No cameras'' data available.\n');
    return;
end;

ind_active_views = find(active_imgviews(:)');
% The unit of the drawing reference: nu
nu = max(lamda);

if ~exist('show_camera','var'),
    show_camera = 1;
end;

figure(4);
clf;
hold on;
if show_camera,
    BASE = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]*nu;  % 原点 x基矢量 原点 y基矢量 原点 z基矢量
    DeltaXYZ = [BASE(:,[2,4,6])*(1+delta), -[1;1;1]*nu*delta];     % text位置： |x, y, z, o

    % 求出图片像素坐标的四个顶点(0,0), (nx-1,0), (nx-1,ny-1), (0,ny-1)的摄像机归一化坐标
    % 若xp和xd都是齐次坐标(3*N, 第三行为1), 则有 xp= KK*xd; 所以xd=KK\xp.
    % draw camers: camera axes base vector initialization (wrt reference frame)
    % 按照Matlab坐标系特点，将摄像机坐标系基矢量变为(e1, e3, -e2),相当于绕e1旋转-pi/2，使深度方向水平
    active_view = find(any(active_imgviews,2)');
    for pp = active_view,
        fc = fc_mat(:,pp);
        cc = cc_mat(:,pp);
        alpha_c = alpha_vec(pp);
        nx = imsize(1,pp);
        ny = imsize(2,pp);
        KK = [fc(1) fc(1)*alpha_c cc(1);0 fc(2) cc(2); 0 0 1];
        IP = KK\[0 nx-1 nx-1 0 0; 0 0 ny-1 ny-1 0; 1 1 1 1 1]*nu;
        IP = reshape([IP;BASE(:,1)*ones(1,5);IP],3,15);
        % Change of reference: wrt reference frame
        Rckk = rodrigues(Omcw(:,pp));
        if hand_list(pp)~=1,
            Rckk(:,3)=-Rckk(:,3);
        end;
        Twkk = Tcw(:,pp);
        BASE_kk = Rckk'*(BASE-Twkk(:,ones(1,6)));       % Xw=Rk'*(Xk-Tk)
        DELTA_kk = Rckk'*(DeltaXYZ-Twkk(:,ones(1,4)));
        IP_kk = Rckk'*(IP -Twkk(:,ones(1,15)));
        figure(4);
        plot3(BASE_kk(1,:),BASE_kk(3,:),-BASE_kk(2,:),'b-','linewidth',1);
        plot3(IP_kk(1,:),IP_kk(3,:),-IP_kk(2,:),'r-','linewidth',1);
        text(DELTA_kk(1,1),DELTA_kk(3,1),-DELTA_kk(2,1),'X','HorizontalAlignment','center');
        text(DELTA_kk(1,2),DELTA_kk(3,2),-DELTA_kk(2,2),'Y','HorizontalAlignment','center');
        text(DELTA_kk(1,3),DELTA_kk(3,3),-DELTA_kk(2,3),'Z','HorizontalAlignment','center');
        text(DELTA_kk(1,4),DELTA_kk(3,4),-DELTA_kk(2,4),['Cam-' num2str(pp)],'HorizontalAlignment','center');
    end;
end;

% draw 3D pts:
step = 100; % step to note the stick
if exist('Xrod','var')==1,
    xyz = permute(reshape(Xrod,3,np1D,n_ima), [2,3,1]);
    plot3(xyz(:,:,1),xyz(:,:,3),-xyz(:,:,2), '.-');
    for i=1:step:n_ima,
        text(Xori(1,i),Xori(3,i),-Xori(2,i),num2str(i),'fontsize',10,'color','k','horizontalalignment','center');
    end;
end;


az = 50;
el = 20;
figure(4);
rotate3d on;
grid on;
title('Extrinsic results of the camera system');
set(4,'color',[1 1 1],'Name','3D','NumberTitle','off');
axis equal vis3d tight;
xlabel('X');
ylabel('Z');
zlabel('-Y');
view(az,el);
axis([Xmin,Xmax,Zmin,Zmax,-Ymax,-Ymin]*1.2);
hold off;
