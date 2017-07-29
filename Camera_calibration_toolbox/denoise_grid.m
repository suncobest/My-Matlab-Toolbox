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
% �Ե��ƽ�������ֵ�ֽ⣬����ֵ��С��������Ӧ����ʸ�������ϵ�����̶ȡ�
% �������ľࣺ�����Ƽ�ȥƽ��ֵ����ƽ�Ƶ����Ƽ������Ĵ�������Э�������
% ����ԭ��أ�ֱ�ӶԵ����������

% ��С���˷�����ֱ��lx_k
for k=1:N_x,
    X = [x_grid(k,:);y_grid(k,:);ones(1,N_y)];
    X = X*X';
    [~,~,V] = svd(X);  % ��ƽ������ȷ���ƽ��1�����࿴��1��ԭ��ƽ���ԭ����Ϊ����
    vert(:,k) = V(:,3);
end;
% X=[x_grid(k,:);y_grid(k,:);ones(1,N_y)]��ʾ��k������ǵ����꣬��ֱ��lx_k�ϡ�
% ��С����ֵ��Ӧ������ʸ��Ϊֱ��lx_k��ͶӰ��������ƽ��ķ�ʸ�����˷����Ӧ��ֱ�ߵ�ת���������
% vert��ʾN_x��ֱ������������ƽ��ķ���������Щ������ȫ����ֱ����ʧ��vx�����ĵ����ߣ������ĵ�vx�ĵ�λʸ��ΪV_vert��

% ��С���˷�����ֱ��ly_k
for k=1:N_y,
    X = [x_grid(:,k)';y_grid(:,k)';ones(1,N_x)];
    X = X*X';
    [~,~,V] = svd(X);
    hori(:,k) = V(:,3);
end;
% X=[x_grid(:,k);y_grid(:,k);ones(1,N_x)]��ʾ��k������ǵ����꣬��ֱ��ly_k��
% hori��ʾN_y��ֱ������������ƽ��ķ���������Щ������ȫ����ֱ����ʧ��vy�����ĵ����ߣ������ĵ�vy�ĵ�λʸ��ΪV_hori��

% 2 principle Vanishing points:������ָ����ʧ��ĵ�λʸ����
X = vert*vert';
[~,~,V] = svd(X);
V_vert = V(:,3);                % ��С���˷�����ʧ�㣨�ǵ���������
X = hori*hori';
[~,~,V] = svd(X);
V_hori = V(:,3);
% �����ģ���ƽ��ƽ������ά���������ƽ�棬��ƽ������ƽ��Ľ���Ϊ��ʧ�ߣ�X��Y�������ʧ�㶼����ʧ���ϡ�
% ���е�lx������ʧ��vx�����е�ly������ʧ��vy��vx�����ĵ����������й����ĺ�ֱ��lx��ƽ���ϣ���˴�ֱ������ƽ��ķ���
% ������ƽ��ķ�����vert��������ֵ�ֽ⣬����Щ�������������Ĳ���ֱ��V_vert������С����ֵ��Ӧ������ʸ������V_vert

% ��ƽ��Ľ��ߴ�ֱ����ƽ�淨ʸ������px_k��ʾ�����ĺ�ֱ��lx_k��ƽ�棬py_k��ʾ�����ĺ�ֱ��ly_k��ƽ��
% ƽ��Px_k�ķ�����Ϊֱ��lx_k�������꣬ƽ��py_k�ķ�����Ϊֱ��ly_k��������
% ��ƽ��Ľ��߶�Ӧ��ֱ�ߵĽ���

for k=1:N_x,
    lx_k = vert(:,k) - dot(V_vert,vert(:,k)) * V_vert;        % ���������������Ӱ�졣��������£�V_vert��ֱ��vert�������С�(��֤���е�lxֱ�߹���ʧ��V_vert)
    vert(:,k) = lx_k;
end;
for k=1:N_y,
    ly_k = hori(:,k) - dot(V_hori,hori(:,k)) * V_hori;        % ���������������Ӱ�졣��������£�V_hori��ֱ��hori�������С�(��֤���е�lyֱ�߹���ʧ��V_hori)
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

x1 = ideal_grid(:,1);                 % ƽ��py_1��ƽ��px_1�Ľ��߷��򣬴�����ָ���һ���ǵ�
x2 = ideal_grid(:,N_x);               % ƽ��py_1��ƽ��px_N_x�Ľ��߷��򣬴�����ָ��x�������һ���ǵ�
x3 = ideal_grid(:,end);               % ƽ��py_N_y��ƽ��px_N_x�Ľ��߷��򣬴�����ָ��ԽǴ���Զ�Ľǵ�
x4 = ideal_grid(:,N_x*(N_y-1)+1);     % ƽ��py_N_y��ƽ��px_1�Ľ��߷��򣬴�����ָ��y�������һ���ǵ�


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

