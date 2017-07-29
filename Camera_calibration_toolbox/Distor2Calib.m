function [fc_2,Rc_2,Tc_2,H_2,distance,V_vert,V_hori,x_all_c,V_hori_pix,V_vert_pix,V_diag1_pix,V_diag2_pix]=Distor2Calib(k_dist,grid_pts_centered,n_sq_x,n_sq_y,Np,W,L,Xgrid_2,f_ini,N_iter,two_focal)

% Computes the calibration parameters knowing the
% distortion factor k_dist

% grid_pts_centered are the grid point coordinates after substraction of
% the optical center.

% can give an optional guess for the focal length f_ini (can set to [])
% can provide the number of iterations for the Iterative Vanishing Point Algorithm

% if the focal length is known perfectly, then, there is no need to iterate,
% and therefore, one can fix: N_iter = 0;

% California Institute of Technology
% (c) Jean-Yves Bouguet - October 7th, 1997
% Beijing Institute of Technology
% Zhang Pengfei - May 4th, 2015


if nargin < 11,
    two_focal = 0;
    if nargin < 10,
        N_iter = 10;
        if nargin < 9,
            f_ini = [1;1];
        end;
    end;
end;

if length(f_ini)<2,
    f_ini=[f_ini;f_ini];
end;

N_x = n_sq_x+1;
N_y = n_sq_y+1;

x_grid = zeros(N_x,N_y);
y_grid = zeros(N_x,N_y);


%%% Computation of the four vanishing points in pixels
x_grid(:) = grid_pts_centered(1,:);   % 以像平面中心为原点的坐标系
y_grid(:) = grid_pts_centered(2,:);

vert = zeros(3,N_x);
hori = zeros(3,N_y);

% [u,s,v] = svd(X), X=u*s*v'; X*X'= u*s^2*u'
% X*X' autocorrelation
% 对点云进行奇异值分解，奇异值大小代表在相应奇异矢量方向上的延伸程度。
% 二阶中心距：将点云减去平均值，即平移到点云几何中心处，再求协方差矩阵。
% 二阶原点矩：直接对点云求自相关

for k=1:N_x,
    X = [x_grid(k,:);y_grid(k,:);ones(1,N_y)];
    X = X*X';
    [~,~,V] = svd(X);    % 像平面向深度方向平移1，焦距看做1，原像平面中心为原点，视为射心（初始位置）
    vert(:,k) = V(:,3);
end;
% X=[x_grid(k,:);y_grid(k,:);ones(1,N_y)]表示第k行网格角点坐标，在直线lx_k上。
% 最小奇异值对应的奇异矢量为直线lx_k和投影中心所在平面的法矢量（此方向对应的直线点转动惯量最大）
% vert表示N_x条直线与射心所在平面的法向量，这些法向量全部垂直于消失点vx和射心的连线（从射心到vx的单位矢量为V_vert）


for k=1:N_y,
    X = [x_grid(:,k)';y_grid(:,k)';ones(1,N_x)];
    X = X*X';
    [~,~,V] = svd(X);
    hori(:,k) = V(:,3);
end;
% X=[x_grid(:,k);y_grid(:,k);ones(1,N_x)]表示第k列网格角点坐标，在直线ly_k上
% hori表示N_y条直线与射心所在平面的法向量，这些法向量全部垂直于消失点vy和射心的连线（从射心到vy的单位矢量为V_hori）

% 2 principle Vanishing points:（主点指向消失点的单位矢量）
X = vert*vert';
[~,~,V] = svd(X);
V_vert = V(:,3);                % 最小二乘法求消失点（角点有噪声）
X = hori*hori';
[~,~,V] = svd(X);
V_hori = V(:,3);
% 过射心，作平面平行于三维世界的网格平面，新平面与像平面的交线为消失线，X和Y方向的消失点都在消失线上。
% 所有的lx交于消失点vx，所有的ly交于消失点vy。vx和射心的连线在所有过射心和直线lx的平面上，因此垂直于所有平面的法线
% 对所有平面的法向量vert进行奇异值分解，且这些法向量都过射心并垂直于V_vert，则最小奇异值对应的奇异矢量就是V_vert

% 两平面的交线垂直于两平面法矢量，令px_k表示过射心和直线lx_k的平面，py_k表示过射心和直线ly_k的平面
% 平面Px_k的法向量为直线lx_k的线坐标，平面py_k的法向量为直线ly_k的线坐标
% 两平面的交线对应两直线的交点

% Square warping:
vert_first = vert(:,1) - dot(V_vert,vert(:,1)) * V_vert;        % 正交化以消除误差影响。理想情况下，V_vert垂直于vert的所有列。第一列为lx_1和射心所在平面的法矢量
vert_last = vert(:,N_x) - dot(V_vert,vert(:,N_x)) * V_vert;     % 最后一列为lx_N_x和射心所在平面的法矢量

hori_first = hori(:,1) - dot(V_hori,hori(:,1)) * V_hori;        % 正交化以消除误差影响。理想情况下，V_hori垂直于hori的所有列。第一列为ly_1和射心所在平面的法矢量
hori_last = hori(:,N_y) - dot(V_hori,hori(:,N_y)) * V_hori;     % 最后一列为ly_N_y和射心所在平面的法矢量


x1 = cross(hori_first,vert_first);   % 平面py_1和平面px_1的交线方向，从射心指向第一个角点
x2 = cross(hori_first,vert_last);    % 平面py_1和平面px_N_x的交线方向，从射心指向x轴上最后一个角点
x3 = cross(hori_last,vert_last);     % 平面py_N_y和平面px_N_x的交线方向，从射心指向对角处最远的角点
x4 = cross(hori_last,vert_first);    % 平面py_N_y和平面px_1的交线方向，从射心指向y轴上最后一个角点

% 将各矢量化为归一化平面上的点，此时焦距为1，因此各点在像平面上，这相当于对四个角点的优化（同样的方法可以对所有角点优化，优化后结果严格满足几何约束）
x1 = x1/x1(3);
x2 = x2/x2(3);
x3 = x3/x3(3);
x4 = x4/x4(3);

square = Rectangle2Square([x1 x2 x3 x4],W,L);     % 将恢复出来的世界平面上的网格变成正方形网格（x1为网格原点），并投影到归一化平面上,得到正方形网格的像

y1 = square(:,1);
y2 = square(:,2);
y3 = square(:,3);
y4 = square(:,4);

H2 = cross(V_vert,V_hori);   % 世界网格平面的法线方向，垂直于dx与dy的方向，也是世界网格平面的消失线

V_diag1 = cross(cross(y1,y3),H2);   % 对角线上的消失点
V_diag2 = cross(cross(y2,y4),H2);

V_diag1 = V_diag1 / norm(V_diag1);  % 从射心指向对角线消失点的单位矢量
V_diag2 = V_diag2 / norm(V_diag2);

V_hori_pix = V_hori;
V_vert_pix = V_vert;
V_diag1_pix = V_diag1;
V_diag2_pix = V_diag2;


% end of computation of the vanishing points in pixels.


fc_2 = f_ini;
x_all_c = [grid_pts_centered(1,:)/fc_2(1);grid_pts_centered(2,:)/fc_2(2)];
x_all_c = comp_distortion(x_all_c,k_dist); % we can this time!!!

% The vertical lines: vert, Horizontal lines: hori
vert = zeros(3,N_x);
hori = zeros(3,N_y);

% 直线和点满足方程x'*l=0，齐次方程的解为矩阵x'最小奇异值对应的右奇异矢量，即矩阵x的左奇异矢量
for counter_k = 1:N_iter, 	% the Iterative Vanishing Points Algorithm to
    % estimate the focal length accurately
    
    x_grid(:) = x_all_c(1,:);
    y_grid(:) = x_all_c(2,:);
    
    for k=1:N_x,
        X = [x_grid(k,:);y_grid(k,:);ones(1,N_y)];
        X = X*X';
        [~,~,V] = svd(X);
        vert(:,k) = V(:,3);
    end;
    
    for k=1:N_y,
        X = [x_grid(:,k)';y_grid(:,k)';ones(1,N_x)];
        X = X*X';
        [~,~,V] = svd(X);
        hori(:,k) = V(:,3);
    end;
    
    % 2 principle Vanishing points:
    X = vert*vert';
    [~,~,V] = svd(X);
    V_vert = V(:,3);
    X = hori*hori';
    [~,~,V] = svd(X);
    V_hori = V(:,3);
   
    
    
    % Square warping:  
    vert_first = vert(:,1) - dot(V_vert,vert(:,1)) * V_vert;
    vert_last = vert(:,N_x) - dot(V_vert,vert(:,N_x)) * V_vert;
    
    hori_first = hori(:,1) - dot(V_hori,hori(:,1)) * V_hori;
    hori_last = hori(:,N_y) - dot(V_hori,hori(:,N_y)) * V_hori;
    
    
    x1 = cross(hori_first,vert_first);
    x2 = cross(hori_first,vert_last);
    x3 = cross(hori_last,vert_last);
    x4 = cross(hori_last,vert_first);
    
    x1 = x1/x1(3);
    x2 = x2/x2(3);
    x3 = x3/x3(3);
    x4 = x4/x4(3);
    
    
    
    [square] = Rectangle2Square([x1 x2 x3 x4],W,L);
    
    y1 = square(:,1);
    y2 = square(:,2);
    y3 = square(:,3);
    y4 = square(:,4);
    
    H2 = cross(V_vert,V_hori);
    
    V_diag1 = cross(cross(y1,y3),H2);
    V_diag2 = cross(cross(y2,y4),H2);
    
    V_diag1 = V_diag1 / norm(V_diag1);
    V_diag2 = V_diag2 / norm(V_diag2);
    
    
    % Estimation of the focal length, and normalization:
    
    % Compute the ellipsis of (1/f^2) positions:
    % a * (1/fx)^2 + b * (1/fx)^2 = -c
    
    
    a1 = V_hori(1);
    b1 = V_hori(2);
    c1 = V_hori(3);
    
    a2 = V_vert(1);
    b2 = V_vert(2);
    c2 = V_vert(3);
    
    a3 = V_diag1(1);
    b3 = V_diag1(2);
    c3 = V_diag1(3);
    
    a4 = V_diag2(1);
    b4 = V_diag2(2);
    c4 = V_diag2(3);
    
    
    if two_focal,
        A = [a1*a2 b1*b2;a3*a4 b3*b4];
        b = -[c1*c2;c3*c4];
        f = sqrt(abs(1./(inv(A)*b)));
    else
        f = sqrt(abs(-(c1*c2*(a1*a2 + b1*b2) + c3*c4*(a3*a4 + b3*b4))/(c1^2*c2^2 + c3^2*c4^2)));
        f = [f;f];
    end;
    
    
    
    % REMARK:
    % if both a and b are small, the calibration is impossible.
    % if one of them is small, only the other focal length is observable
    % if none is small, both focals are observable
    fc_2 = fc_2 .* f;
    
    
    % DEBUG PART: fix focal to 500...
    %fc_2= [500;500]; disp('Line 293 to be earased in Distor2Calib.m');
    
    
    % end of focal compensation
    
    % normalize by the current focal:
    
    x_all = [grid_pts_centered(1,:)/fc_2(1);grid_pts_centered(2,:)/fc_2(2)];
    
    % Compensate by the distortion factor:
    
    x_all_c = comp_distortion(x_all,k_dist);
    
end;

% At that point, we hope that the distortion is gone...

x_grid(:) = x_all_c(1,:);
y_grid(:) = x_all_c(2,:);

for k=1:N_x,
    X = [x_grid(k,:);y_grid(k,:);ones(1,N_y)];
    X = X*X';
    [~,~,V] = svd(X);
    vert(:,k) = V(:,3);
end;

for k=1:N_y,
    X = [x_grid(:,k)';y_grid(:,k)';ones(1,N_x)];
    X = X*X';
    [~,~,V] = svd(X);
    hori(:,k) = V(:,3);
end;

% Vanishing points:
    X = vert*vert';
    [~,~,V] = svd(X);
    V_vert = V(:,3);
    X = hori*hori';
    [~,~,V] = svd(X);
    V_hori = V(:,3);

% Horizon:

H_2 = cross(V_vert,V_hori);


% pick a plane in front of the camera (positive depth)
if H_2(3) < 0, H_2 = -H_2; end;


% Rotation matrix:

if V_hori(1) < 0, V_hori = -V_hori; end;

V_hori = V_hori/norm(V_hori);
H_2 = H_2/norm(H_2);

V_hori = V_hori - dot(V_hori,H_2)*H_2;

Rc_2 = [V_hori cross(H_2,V_hori) H_2];

Rc_2 = Rc_2 / det(Rc_2);

%omc_2 = rodrigues(Rc_2);

%Rc_2 = rodrigues(omc_2);

% Find the distance of the plane for translation vector:

xc_2 = [x_all_c;ones(1,Np)];

Zc_2 = 1./sum(xc_2 .* (Rc_2(:,3)*ones(1,Np)));

Xo_2 = [sum(xc_2 .* (Rc_2(:,1)*ones(1,Np))).*Zc_2 ; sum(xc_2 .* (Rc_2(:,2)*ones(1,Np))).*Zc_2];

XXo_2 = Xo_2 - mean(Xo_2')'*ones(1,Np);

distance_x = norm(Xgrid_2(1,:))/norm(XXo_2(1,:));
distance_y = norm(Xgrid_2(2,:))/norm(XXo_2(2,:));


distance = sum(sum(XXo_2(1:2,:).*Xgrid_2(1:2,:)))/sum(sum(XXo_2(1:2,:).^2));

alpha = abs(distance_x - distance_y)/distance;

if (alpha>0.1)&&~two_focal,
    disp('Should use two focals in x and y...');
end;

% Deduce the translation vector:

Tc_2 = distance * H_2;


return;


