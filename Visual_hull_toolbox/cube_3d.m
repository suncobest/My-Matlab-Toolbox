d_sq = 30; % cube edge length
cube = [0 1 1 0 0 1 1 0
        0 0 1 1 0 0 1 1
        0 0 0 0 1 1 1 1] * d_sq;    % cube vertices 

Xw = cube;
% Xwl = [Xw(:,1:4),Xw(:,1),Xw(:,5:8),Xw(:,5:6),Xw(:,2:3),Xw(:,7:8),Xw(:,4)]; % cube edge

n_step = 3;  % number of motion steps

om1 = [0;0;0];
T1 = -d_sq/2 * [1;1;1];

om2 = [2;-1;3 ]; % randn(3,1);
T2 = [0;0;0];

om3 = [0;0;0];
T3 = -(d_sq/2+20) * [1;-1;0]+[0;0;930];

om_motion = om1;
T_motion = T1;

%%% compose motion

% for kk = 1:n_step
%     eval(['Xw = rigid_motion(Xw,om' num2str(kk) ',T' num2str(kk) ');']);
% end

for kk = 2:n_step
    eval(['[om_motion,T_motion] = compose_motion(om_motion,T_motion,om' num2str(kk) ',T' num2str(kk) ');']);
end

Xw = rigid_motion(Xw,om_motion,T_motion);

save Xw_3d.mat d_sq cube Xw om_motion T_motion 

clear


%% 以Cam1的坐标系为世界坐标系，对于已知的三维点Xw，求出其在Camkk下的摄像机坐标Xckk，和投影的像素坐标xpkk

load('Xw_3d.mat');

path_dir = 'D:\My Documents\MATLAB\Camera_Calibration\Calibration\131008A\';
load([path_dir,'Camera_list.mat'])
N = size(Xw,2);

for kk = 1:n_cam
    
    % Relative position of camera_kk wrt camera_1: (om,T,hand)
    R = rodrigues(om_list(:,kk)) * diag([1 1 hand_list(kk)]);
    eval(['Xc' num2str(kk) ' = R * Xw + repmat(T_list(:,kk),[1 N]);']);
    eval(['xp' num2str(kk) ' = project_points_mirror(Xw,om_list(:,kk),T_list(:,kk),hand_list(kk),fc_list(:,kk),cc_list(:,kk),kc_list(:,kk),alpha_c_list(:,kk));']);
end

%% 由每个摄像机的像素坐标xpkk和摄像机参数Camera_list.mat求出在Camkk下的摄像机坐标Xckk
% 若xp和xd都是齐次坐标（3*N，第三行为1），则有
% xp = [1 0 cc(1);0 1 cc(2);0 0 1]*[fc(1) 0 0;0 fc(2) 0;0 0 1]*[1 alpha_c 0;0 1 0;0 0 1] * xd = KK * xd
% 其中KK = [fc(1) alpha_c*fc(1) cc(1);0 fc(2) cc(2);0 0 1]

% N = size(xp1,2);
fc_left = fc_list(:,1);
cc_left = cc_list(:,1);
kc_left = kc_list(:,1);
alpha_c_left = alpha_c_list(:,1);

XL_list = zeros(3*N,n_cam-1);
for kk = 2:n_cam
    eval(['xR = xp' num2str(kk) ';']);
    om = om_list(:,kk);
    T = T_list(:,kk);
    hand = hand_list(kk);
    
    fc_right=fc_list(:,kk);
    cc_right=cc_list(:,kk);
    kc_right=kc_list(:,kk);
    alpha_c_right=alpha_c_list(:,kk);
    
    [XL,XR] = stereo_triangulation2(xp1,xR,om,T,hand,fc_left,cc_left,kc_left,alpha_c_left,fc_right,cc_right,kc_right,alpha_c_right);
    XL_list(:,kk-1) = XL;   % Xc被重复计算了n_cam-1遍，所以XL_list有n_cam-1列，3*N行
    eval(['Yc' num2str(kk) '= XR;']);
end

Yc1 =  mean(XL_list,2);  % 对n_cam-1个Xc取平均（沿列标增长的方向取平均）
Yc1_std = std(XL_list,0,2);  % take the normalization by N-1(default)
sigma3 = abs(XL_list - Yc1(:,ones(1,n_cam-1))) > 3 * Yc1_std(:,ones(1,n_cam-1));
if any(any(sigma3))
    warning('Error is too large!');
end
Yc1 = reshape(Yc1,3,N);
