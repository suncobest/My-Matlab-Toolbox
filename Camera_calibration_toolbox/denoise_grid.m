function [ideal_grid,V_hori,V_vert,H]=denoise_grid(grid_pts,n_sq_x,n_sq_y)

% Computes the ideal grid and vanishing points giving the corner points possible with noise. 
%
% grid_pts are the grid points with noise.
% n_sq_x is number of square in the 1st dimension
% n_sq_y is number of square in the 2nd dimension
%
% ideal_grid are the denoise grid points.
% V_hori is the vanishing points along the 1st dimension.
% V_vert is the vanishing points along the 2nd dimension.
% H is the homography which map a unit square to the four outer points of grid_pts.


N_x = n_sq_x+1;
N_y = n_sq_y+1;

[mx,Np] = size(grid_pts); 
assert(mx==2||mx==3,'Unexpected points dimmension!');
assert(Np==N_x*N_y,'Number of points unmatched!');

if mx == 3,
    grid_pts = grid_pts./(ones(3,1)*grid_pts(3,:));
end;


x_grid = zeros(N_x,N_y);
y_grid = zeros(N_x,N_y);
x_grid(:) = grid_pts(1,:); 
y_grid(:) = grid_pts(2,:);

vert = zeros(3,N_x);
hori = zeros(3,N_y);

% [u,s,v] = svd(X), X=u*s*v'; X*X'= u*s^2*u'
% X*X' autocorrelation
% 对点云进行奇异值分解，奇异值大小代表在相应奇异矢量方向上的延伸程度。
% 二阶中心距：将点云减去平均值，即平移到点云几何中心处，再求协方差矩阵。
% 二阶原点矩：直接对点云求自相关

% 最小二乘法计算直线lx_k
for k=1:N_x,
    X = [x_grid(k,:);y_grid(k,:);ones(1,N_y)];
    X = X*X';
    [~,~,V] = svd(X);  % 像平面向深度方向平移1，焦距看做1，原像平面的原点作为射心
    vert(:,k) = V(:,3);
end;
% X=[x_grid(k,:);y_grid(k,:);ones(1,N_y)]表示第k行网格角点坐标，在直线lx_k上。
% 最小奇异值对应的奇异矢量为直线lx_k和投影中心所在平面的法矢量（此方向对应的直线点转动惯量最大）
% vert表示N_x条直线与射心所在平面的法向量，这些法向量全部垂直于消失点vx和射心的连线（从射心到vx的单位矢量为V_vert）

% 最小二乘法计算直线ly_k
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

for k=1:N_x,
    lx_k = vert(:,k) - dot(V_vert,vert(:,k)) * V_vert;        % 正交化以消除误差影响。理想情况下，V_vert垂直于vert的所有列。(保证所有的lx直线过消失点V_vert)
    vert(:,k) = lx_k;
end;
for k=1:N_y,
    ly_k = hori(:,k) - dot(V_hori,hori(:,k)) * V_hori;        % 正交化以消除误差影响。理想情况下，V_hori垂直于hori的所有列。(保证所有的ly直线过消失点V_hori)
    hori(:,k) = ly_k;
end;

% regenerate grid points without error
ideal_grid = zeros(2,Np);    
for j = 1:N_y,
    for i=1:N_x,
        x_ij = cross(hori(:,j),vert(:,i));   
        x_ij = x_ij(1:2)/x_ij(3);
        ideal_grid(:,(j-1)*N_x+i) = x_ij;
    end;
end;

x1 = ideal_grid(:,1);                 % 平面py_1和平面px_1的交线方向，从射心指向第一个角点
x2 = ideal_grid(:,N_x);               % 平面py_1和平面px_N_x的交线方向，从射心指向x轴上最后一个角点
x3 = ideal_grid(:,end);               % 平面py_N_y和平面px_N_x的交线方向，从射心指向对角处最远的角点
x4 = ideal_grid(:,N_x*(N_y-1)+1);     % 平面py_N_y和平面px_1的交线方向，从射心指向y轴上最后一个角点


if abs(V_vert(3))>1e-5,
    V_vert = V_vert/ V_vert(3);
end;

if abs(V_hori(3))>1e-5,
    V_hori = V_hori/ V_hori(3);
end;

Square = [0 1 0 1;0 0 1 1];
H = compute_homography_lm([x1 x2 x4 x3],Square);  % the homography to map a unit square to the image outer grid.


return;






%% test
Nx = 5;
Ny = 5;
scale = randi(100);
[Ygrid,Xgrid]=meshgrid(((1:Ny)-1)*scale,((1:Nx)-1)*scale);
X0=[Xgrid(:)';Ygrid(:)';ones(1,Nx*Ny)];
K=[scale 0 scale/2;0 scale scale/2;0 0 1];
T=rand(3,1)*10*scale;
T=T/sign(T(3));
R=rodrigues(rand(3,1)*0.6);
H=K*[R(:,1:2),T];
x=H*X0;
x=x./(ones(3,1)*x(3,:));
x1234=x(:,[1,Nx,Nx*Ny,Nx*(Ny-1)+1]);
[du,dv]=UnWarpPlane(x1234);
dot(du,dv)

y=x(1:2,:)+rand(2,Nx*Ny)/10;
[ideal_grid,V_hori,V_vert,H]=denoise_grid(x,Nx-1,Ny-1);
xgrid = zeros(Nx,Ny);
ygrid = zeros(Nx,Ny);
xgrid(:) = ideal_grid(1,:);
ygrid(:) = ideal_grid(2,:);

figure(1);
plot(x(1,:),x(2,:),'ro')            % original data
hold on,axis equal,plot(y(1,:),y(2,:),'b+')   % with noise
plot(ideal_grid(1,:),ideal_grid(2,:),'k.')    % denoise

if V_vert(3)==1 && V_hori(3)==1 && max(abs([V_vert;V_hori]))<10*max(abs(ideal_grid(:))),
    xhori_line = [xgrid([1,Nx],:);ones(1,Ny)*V_hori(1)];
    yhori_line = [ygrid([1,Nx],:);ones(1,Ny)*V_hori(2)];
    xvert_line = [xgrid(:,[1,Ny])';ones(1,Nx)*V_vert(1)];
    yvert_line = [ygrid(:,[1,Ny])';ones(1,Nx)*V_vert(2)];
    plot(xhori_line,yhori_line,'r-',V_hori(1),V_hori(2),'ro')       % plot(xhori_line,yhori_line,'r-')
    plot(xvert_line,yvert_line,'r-',V_vert(1),V_vert(2),'ro')       % plot(xvert_line,yvert_line,'r-')
else
    plot(xgrid,ygrid,'r-')
    plot(xgrid',ygrid','r-')
end

[ideal_grid,V_hori,V_vert,H]=denoise_grid(y,Nx-1,Ny-1);
xgrid = zeros(Nx,Ny);
ygrid = zeros(Nx,Ny);
xgrid(:) = ideal_grid(1,:);
ygrid(:) = ideal_grid(2,:);
plot(ideal_grid(1,:),ideal_grid(2,:),'go')    % denoise
if V_vert(3)==1 && V_hori(3)==1 && max(abs([V_vert;V_hori]))<10*max(abs(ideal_grid(:))),
    xhori_line = [xgrid([1,Nx],:);ones(1,Ny)*V_hori(1)];
    yhori_line = [ygrid([1,Nx],:);ones(1,Ny)*V_hori(2)];
    xvert_line = [xgrid(:,[1,Ny])';ones(1,Nx)*V_vert(1)];
    yvert_line = [ygrid(:,[1,Ny])';ones(1,Nx)*V_vert(2)];
    plot(xhori_line,yhori_line,'g-',V_hori(1),V_hori(2),'go')       % plot(xhori_line,yhori_line,'g-')
    plot(xvert_line,yvert_line,'g-',V_vert(1),V_vert(2),'go')       % plot(xvert_line,yvert_line,'g-')
else
    plot(xgrid,ygrid,'g-')
    plot(xgrid',ygrid','g-')
end

axis tight

