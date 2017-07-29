function [Qckk,Tckk] = compute_extrinsic_init_mirror(x_kk,X_kk,fc,cc,kc,alpha_c,hand)

% compute_extrinsic_init_mirror
%
% Computes the extrinsic parameters attached to a 3D structure X_kk given its projection
% on the image plane x_kk and the intrinsic camera parameters fc, cc and kc.
% Works with planar and non-planar structures.
%
% INPUT: x_kk: Feature locations on the images
%       X_kk: Corresponding grid coordinates
%       fc: Camera focal length
%       cc: Principal point coordinates
%       kc: Distortion coefficients
%       alpha_c: Skew coefficient
%       hand: Handness of camera frame wrt world frame (1 or -1);
%
% OUTPUT: Qckk: 4D quaternion rotation vector attached to the grid positions in space
%              Tckk: 3D translation vector attached to the grid positions in space
%
% Method: Computes the normalized point coordinates, then computes the 3D pose
%
% Important functions called within that program:
% normalize_pixel: Computes the normalize image point coordinates.
% compute_homography_lm: Computes homography with lm iteration algorithm.
% trans_quat_mat: transform rotation matrix to quternion vector, or vice versa.
%
% See also compute_extrinsic_init.

if nargin < 7,
    hand = 1;
    if nargin < 6,
        alpha_c = 0;
        if nargin < 5,
            kc = zeros(5,1);
            if nargin < 4,
                cc = zeros(2,1);
                if nargin < 3,
                    fc = ones(2,1);
                    if nargin < 2,
                        error('Need 2D projections and 3D points (in compute_extrinsic.m)');
                    end;
                end;
            end;
        end;
    end;
end;

%keyboard;
% Compute the normalized coordinates:
xn = normalize_pixel(x_kk,fc,cc,kc,alpha_c);
Np = size(xn,2);

if size(X_kk,1)== 2,
    H = compute_homography_lm(xn,X_kk);  % 计算规一坐标与随体坐标的单应性矩阵   
    % H=[h1 h2 h3]=aK[r1 r2 t],a为待定值，因为xn=[Xc/Zc;Yc/Zc;1],所以K=eye(3),H=[h1 h2 h3]=[r1 r2 t]/Zc
    % De-embed the motion parameters from the homography: 
    sc = mean([norm(H(:,1));norm(H(:,2))]);  % 求出h1，h2的平均长度
    H = H/sc; % 对H进行归一化
    
    if H(9)<0,              % Tc3>0
        H=-H;
    end;
    
    u1 = H(:,1);
    u1 = u1 / norm(u1);
    u2 = H(:,2) - dot(u1,H(:,2)) * u1;
    u2 = u2 / norm(u2);
    u3 = cross(u1,u2);
    Rckk = [u1 u2 u3];
    Tckk = H(:,3);
    Qckk = trans_quat_mat(Rckk);
    return;
end;

% Check for planarity of the structure:
X_mean = mean(X_kk,2);  % 均值, 重心位置
% 误差函数: 各点相对重心的位置，X_mean*ones(1,Np)=X_mean(:,ones(1,Np))
Y = X_kk - (X_mean*ones(1,Np));  
YY = Y*Y'; % 协方差矩阵
% 对于三维情形来说，协方差矩阵与转动惯量密切相关。若Np个点的总质量m=1,
% 则转动惯量I=(tr(Y'*Y)*eye(3)-Y*Y')/Np。
% 其中tr(Y'*Y)表示矩阵Y'*Y的迹，等于矩阵Y的所有元素平方和，tr(Y'*Y)=Y(:)'*Y(:)。
% 主成分分析PCA， YY*V=V*S，V的列向量为点系的主轴，
% S=X_new*X_new',其中X_new=V'*(X_kk-X_mean)
[~,S,V] = svd(YY); 
r = S(3,3)/S(2,2);  % 最小奇异值与第二小奇异值的比值

if (r < 1e-3)||(Np < 4), %1e-3, %1e-4, %norm(X_kk(3,:)) < eps, % Test of planarity
    %fprintf(1,'Planar structure detected: r=%f\n',r);  
    % Transform the plane to bring it in the Z=0 plane:
    R_transform = V';
    % 酉矩阵V=[v1,v2,v3]，其中v1，v2，v3分别为共面点阵X_kk的三个主方向，且v3为法方向。
    % R_transform=[v1t;v2t;v3t];
    %norm(R_transform(1:2,3))
    
    if norm(R_transform(1:2,3)) < 1e-6,   % check if R_transform(:,3)==[0;0;1]
        R_transform = eye(3);                      % 若R_transform(:,3)==[0;0;1]，则v3=[0;0;1]为z方向
    elseif det(R_transform) < 0, 
        R_transform = -R_transform;         %  YY = V*SV' = -V*S*(-V')
    end;  
    
    T_transform = -(R_transform)*X_mean;
    % X_new为随体坐标（以重心X_mean为原点，主方向v1，v2，v3为基矢量）
    X_new = R_transform*X_kk + T_transform(:, ones(1,Np));  
    
    %  X_kk= V*X_new + X_mean*ones(1,Np); （显然X_new(3,:)=0） 
    % Compute the planar homography:  
    H = compute_homography_lm(xn,X_new(1:2,:));  % 计算规一坐标与随体坐标的单应性矩阵
    sc = mean([norm(H(:,1));norm(H(:,2))]);  % 求出h1，h2的平均长度
    H = H/sc; % 对H进行归一化
    
    if H(9)<0,
        H=-H;
    end;
    
    u1 = H(:,1);
    u1 = u1 / norm(u1);
    u2 = H(:,2) - dot(u1,H(:,2)) * u1;
    u2 = u2 / norm(u2);
    u3 = cross(u1,u2);
    Rckk = [u1 u2 u3];
    Tckk = H(:,3);
    
    % Because X_new = R_transform*X_kk + T_transform， if Xc = Rckk * X_new + Tckk,
    % then Xc = Rckk * R_transform * X_kk + Rckk* T_transform + Tckk    
    Tckk = Tckk + Rckk* T_transform;
    Rckk = Rckk * R_transform;  
    Qckk = trans_quat_mat(Rckk);       
else
    
    %fprintf(1,'Non planar structure detected: r=%f\n',r);
    
    % Computes an initial guess for extrinsic parameters (works for general 3d structure, not planar!!!):
    % The DLT method is applied here!!
    % xn=[r1 r1 r3 t]*[X_kk;1]/Zc
    % 其中X_kk为不共面的三维点，xn为在Zc=1处的像点，两者构成非单应映射,
    % 所以用直接线性变换DLT求系数矩阵[r1 r1 r3 t]/Zc
    if hand~=1,
        X_kk(3,:) = -X_kk(3,:);
    end;
    J = zeros(2*Np,12);
    
    xX = (ones(3,1)*xn(1,:)).*X_kk;
    yX = (ones(3,1)*xn(2,:)).*X_kk;
    
    J(1:2:end,[1 4 7]) = -X_kk';
    J(2:2:end,[2 5 8]) = X_kk';
    J(1:2:end,[3 6 9]) = xX';
    J(2:2:end,[3 6 9]) = -yX';
    J(1:2:end,12) = xn(1,:)';
    J(2:2:end,12) = -xn(2,:)';
    J(1:2:end,10) = -ones(Np,1);
    J(2:2:end,11) = ones(Np,1);
    
    JJ = J'*J;
    [~,~,V] = svd(JJ);
    
    RR = reshape(V(1:9,12),3,3);  % V(1:9,12)=a[r1;r2;r3],a为待定实数
    
    if det(RR) < 0,
        V(:,12) = -V(:,12);
        RR = -RR;
    end;
    
    [Ur,~,Vr] = svd(RR);   
    Rckk = Ur*Vr';
    
    sc = norm(V(1:9,12)) / norm(Rckk(:));
    Tckk = V(10:12,12)/sc;    
    Qckk = trans_quat_mat(Rckk);
    
end;
