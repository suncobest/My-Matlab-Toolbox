function [fc,omc,Tc] = focus_extrinsic_estimate2(grid_pts,Xgrid,cc,f_ini,k_dist,two_focal,N_iter)

% Computes the focal length and the extrinsic parameters knowing the principal point and distortion factor k_dist
%
% grid_pts are the grid pixel coordinates. Xgrid are the world grid points.
% cc is the principal point. Set it to the center of image if you dont't know it.
% fc is the focal length being estimated. omc is the rotation of world plane, Tc is the translation.
%
% First estimate focus without distortion,then run iteration to compensate for distortion.
% The results will be less accurate when distortion is large.
%
% if the focal length is known perfectly, then, there is no need to iterate,
% and therefore, one can fix: N_iter = 0;
%
% See also focus_extrinsic_estimate, center_extrinsic_estimate, Distor2Calib.

% California Institute of Technology
% (c) Jean-Yves Bouguet - October 7th, 1997
% Beijing Institute of Technology
% Zhang Pengfei - May 23th, 2015


if nargin < 7,
    N_iter = 10;
    if nargin < 6,
        two_focal = 0;
        if nargin < 5,
            k_dist = 0;
            if nargin < 4,
                f_ini = [1;1];
            end;
        end;
    end;
end;

[mx,nx] = size(grid_pts);
[mX,nX] = size(Xgrid);
assert(mx==2||mx==3,'Unexpected points dimension!');
assert(mX==2||mX==3,'Unexpected points dimension!');
assert(nx==nX,'Number of 2D and 3D points unmatched!');
assert(length(k_dist)<6 && ~isempty(k_dist),'Unexpected dimension of distortion coefficient!');
Np = nx;

if length(k_dist)>1,
    k = zeros(5,1);
    k(k_dist~=0) = k_dist(k_dist~=0);
    k_dist = k;
    if norm(k(2:5))~=0,
        fprintf(1,'\nNote: High order of distortion coefficient detected! Make sure the distortion coefficient is right!\n');
    end;
end;

if length(f_ini)==1,
    f_ini=[f_ini;f_ini];
end;

if mx == 3,
    grid_pts = grid_pts./(ones(3,1)*grid_pts(3,:));
end;

grid_pts = grid_pts(1:2,:);
grid_pts_centered = grid_pts - cc(:,ones(1,Np));
Xgrid = Xgrid(1:2,:);
fprintf(1,'\nMake sure that the chessboard corner points are regularly arranged, and the first dimension is along the x axis!\n');
fprintf(1,'\nMake sure the first point of the chessboard is the origin of world coordinate!\n');

idx = find(diff(Xgrid(1,:))<0);
N_x = idx(1);
N_y = floor(Np/N_x);
assert(Np==N_x*N_y,'Error: the corner points is badly arranged!');

W = abs(Xgrid(1,N_x)-Xgrid(1,1));               % X方向网格总长度为W, Y方向网格总长度为L
L = abs(Xgrid(2,N_x*(N_y-1)+1)-Xgrid(2,1));

flag = 0;
if norm(k_dist)~=0,
    flag = 1;
end;


fc = f_ini;
x_all_c = [grid_pts_centered(1,:)/fc(1);grid_pts_centered(2,:)/fc(2)];

x_grid = zeros(N_x,N_y);
y_grid = zeros(N_x,N_y);

%%% Computation of the four vanishing points on the normalized plane
x_grid(:) = x_all_c(1,:);   % 以像平面中心为原点的坐标系
y_grid(:) = x_all_c(2,:);

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
    [~,~,V] = svd(X); % 像平面向深度方向平移1，焦距看做1，原像平面中心为原点，视为射心（初始位置）
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

% Estimation of the focal length, and normalization:

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


A = [a1*a2  b1*b2; a3*a4  b3*b4];
b = -[c1*c2;c3*c4];

if two_focal,
    f = sqrt(abs(1./(A\b)));
else
    f = sqrt(abs(1./(sum(A,2)\b))) * ones(2,1);
end;

fc = fc.*f;

% iterate when k_dist~=0
if flag,
    % The vertical lines: vert, Horizontal lines: hori
    vert = zeros(3,N_x);
    hori = zeros(3,N_y);
    for counter_k = 1:N_iter, 	% the Iterative Vanishing Points Algorithm to
        % normalize by the current focal:
        x_all = [grid_pts_centered(1,:)/fc(1);grid_pts_centered(2,:)/fc(2)];

        % Compensate by the distortion factor:
        x_all_c = comp_distortion(x_all,k_dist);
        x_grid(:) = x_all_c(1,:);
        y_grid(:) = x_all_c(2,:);

        % 直线和点满足方程x'*l=0，齐次方程的解为矩阵x'最小奇异值对应的右奇异矢量，即矩阵x的左奇异矢量
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

        H2 = cross(V_vert,V_hori);                  % ideal line of the 3D plane, also the normal
        V_diag1 = cross(cross(y1,y3),H2);
        V_diag2 = cross(cross(y2,y4),H2);

        V_diag1 = V_diag1 / norm(V_diag1);
        V_diag2 = V_diag2 / norm(V_diag2);

        % Estimation of the focal length, and normalization:
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

        A = [a1*a2  b1*b2; a3*a4  b3*b4];
        b = -[c1*c2;c3*c4];

        if two_focal,
            f = sqrt(abs(1./(A\b)));
        else
            f = sqrt(abs(1./(sum(A,2)\b))) * ones(2,1);
        end;

        % REMARK:
        % if both a and b are small, the calibration is impossible.
        % if one of them is small, only the other focal length is observable
        % if none is small, both focals are observable

        fc = fc .* f;
        % end of focal compensation
    end;
end;

% At that point, we hope that the distortion is gone...
x_all_c = [grid_pts_centered(1,:)/fc(1);grid_pts_centered(2,:)/fc(2)];
if flag,
    x_all_c = comp_distortion(x_all_c,k_dist);
end;

x_grid(:) = x_all_c(1,:);
y_grid(:) = x_all_c(2,:);

vert = zeros(3,N_x);
hori = zeros(3,N_y);

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


x1 = cross(hori_first,vert_first);   % 平面py_1和平面px_1的交线方向，从射心指向第一个角点
x2 = cross(hori_first,vert_last);    % 平面py_1和平面px_N_x的交线方向，从射心指向x轴上最后一个角点
x3 = cross(hori_last,vert_last);     % 平面py_N_y和平面px_N_x的交线方向，从射心指向对角处最远的角点
x4 = cross(hori_last,vert_first);    % 平面py_N_y和平面px_1的交线方向，从射心指向y轴上最后一个角点

x1 = x1/x1(3);
x2 = x2/x2(3);
x3 = x3/x3(3);
x4 = x4/x4(3);

H2 = cross(V_vert,V_hori);

% Rotation matrix:
dx = x2-x1;
dx = dx/norm(dx);
dy = x4-x1;
dy = dy/norm(dy);
dz = cross(dx,dy);

if H2(3) < 0, H2 = -H2; end;        % the normal of the 3D plane
H2 = sign(dz(3))*H2;
H2 = H2/norm(H2);        % 消失线的单位坐标。设过射心和消失线的平面为Pov，则H2为Pov的单位法矢量

idx = abs(dx);
idx = find(idx==max(idx));

if V_hori(idx) < 0, V_hori = -V_hori; end;        % the point at the infinity of x axis
V_hori = sign(dx(idx))*V_hori;
V_hori = V_hori - dot(V_hori,H2)*H2;
V_hori = V_hori/norm(V_hori);
Rc = [V_hori, cross(H2,V_hori), H2];
omc = rodrigues(Rc);

% Deduce the translation vector: x1 is the origin
X1 = x1 / dot(H2,x1);    % d1 = dot(H2,x1)为x1到平面Pov的距离， X1 = x1/d1使得X1与Pov的距离为1;
X2 = x2 / dot(H2,x2);    % d2 = dot(H2,x2)为x2到平面Pov的距离， X2 = x2/d2使得X2与Pov的距离为1;
X3 = x3 / dot(H2,x3);    % d3 = dot(H2,x3)为x3到平面Pov的距离， X3 = x3/d3使得X3与Pov的距离为1;
X4 = x4 / dot(H2,x4);    % d4 = dot(H2,x4)为x4到平面Pov的距离， X4 = x4/d4使得X4与Pov的距离为1;

% 显然x1,x2,x3,x4与X1,X2,X3,X4透视对应，且X1,X2,X3,X4所在平面平行于Pov平面（平行于世界平面），且距离为1。

scale = X1(3);     %   scale = x13/d1

X1 = X1/scale;      % 经过缩放后，X1(3)=1,即X1在归一化平面上，若x1,x2,x3,x4都在归一化平面上，则X1与x1重合。
X2 = X2/scale;      % X1,X2,X3,X4所在平面过X1且平行于Pov平面，因此X1,X2,X3,X4相似于世界网格。
X3 = X3/scale;      % 一般四边形x1,x2,x3,x4恢复为平行四边形X1,X2,X3,X4
X4 = X4/scale;

u_hori = X2 - X1;  % 相邻两边为网格矢量
u_vert = X4 - X1;

scale = mean([W/norm(u_hori); L/norm(u_vert)]);
Tc = scale * X1;

return;






% N_x = n_sq_x+1;
% N_y = n_sq_y+1;
% W = abs(Xgrid(1,N_x)-Xgrid(1,1));               % X方向网格总长度为W, Y方向网格总长度为L
% L = abs(Xgrid(2,N_x*(N_y-1)+1)-Xgrid(2,1));
%
% [ideal_grid,V_hori,V_vert,H]=denoise_grid(grid_pts,n_sq_x,n_sq_y);
% x1 = [ideal_grid(:,1);1];                 % 平面py_1和平面px_1的交线方向，从射心指向第一个角点
% x2 = [ideal_grid(:,N_x);1];               % 平面py_1和平面px_N_x的交线方向，从射心指向x轴上最后一个角点
% x3 = [ideal_grid(:,end);1];               % 平面py_N_y和平面px_N_x的交线方向，从射心指向对角处最远的角点
% x4 = [ideal_grid(:,N_x*(N_y-1)+1);1];     % 平面py_N_y和平面px_1的交线方向，从射心指向y轴上最后一个角点
