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

W = abs(Xgrid(1,N_x)-Xgrid(1,1));               % X���������ܳ���ΪW, Y���������ܳ���ΪL
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
x_grid(:) = x_all_c(1,:);   % ����ƽ������Ϊԭ�������ϵ
y_grid(:) = x_all_c(2,:);

vert = zeros(3,N_x);
hori = zeros(3,N_y);

% [u,s,v] = svd(X), X=u*s*v'; X*X'= u*s^2*u'
% X*X' autocorrelation
% �Ե��ƽ�������ֵ�ֽ⣬����ֵ��С��������Ӧ����ʸ�������ϵ�����̶ȡ�
% �������ľࣺ�����Ƽ�ȥƽ��ֵ����ƽ�Ƶ����Ƽ������Ĵ�������Э�������
% ����ԭ��أ�ֱ�ӶԵ����������

for k=1:N_x,
    X = [x_grid(k,:);y_grid(k,:);ones(1,N_y)];
    X = X*X';
    [~,~,V] = svd(X); % ��ƽ������ȷ���ƽ��1�����࿴��1��ԭ��ƽ������Ϊԭ�㣬��Ϊ���ģ���ʼλ�ã�
    vert(:,k) = V(:,3);
end;
% X=[x_grid(k,:);y_grid(k,:);ones(1,N_y)]��ʾ��k������ǵ����꣬��ֱ��lx_k�ϡ�
% ��С����ֵ��Ӧ������ʸ��Ϊֱ��lx_k��ͶӰ��������ƽ��ķ�ʸ�����˷����Ӧ��ֱ�ߵ�ת���������
% vert��ʾN_x��ֱ������������ƽ��ķ���������Щ������ȫ����ֱ����ʧ��vx�����ĵ����ߣ������ĵ�vx�ĵ�λʸ��ΪV_vert��


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

% Square warping:
vert_first = vert(:,1) - dot(V_vert,vert(:,1)) * V_vert;        % ���������������Ӱ�졣��������£�V_vert��ֱ��vert�������С���һ��Ϊlx_1����������ƽ��ķ�ʸ��
vert_last = vert(:,N_x) - dot(V_vert,vert(:,N_x)) * V_vert;     % ���һ��Ϊlx_N_x����������ƽ��ķ�ʸ��

hori_first = hori(:,1) - dot(V_hori,hori(:,1)) * V_hori;        % ���������������Ӱ�졣��������£�V_hori��ֱ��hori�������С���һ��Ϊly_1����������ƽ��ķ�ʸ��
hori_last = hori(:,N_y) - dot(V_hori,hori(:,N_y)) * V_hori;     % ���һ��Ϊly_N_y����������ƽ��ķ�ʸ��


x1 = cross(hori_first,vert_first);   % ƽ��py_1��ƽ��px_1�Ľ��߷��򣬴�����ָ���һ���ǵ�
x2 = cross(hori_first,vert_last);    % ƽ��py_1��ƽ��px_N_x�Ľ��߷��򣬴�����ָ��x�������һ���ǵ�
x3 = cross(hori_last,vert_last);     % ƽ��py_N_y��ƽ��px_N_x�Ľ��߷��򣬴�����ָ��ԽǴ���Զ�Ľǵ�
x4 = cross(hori_last,vert_first);    % ƽ��py_N_y��ƽ��px_1�Ľ��߷��򣬴�����ָ��y�������һ���ǵ�

% ����ʸ����Ϊ��һ��ƽ���ϵĵ㣬��ʱ����Ϊ1����˸�������ƽ���ϣ����൱�ڶ��ĸ��ǵ���Ż���ͬ���ķ������Զ����нǵ��Ż����Ż������ϸ����㼸��Լ����
x1 = x1/x1(3);
x2 = x2/x2(3);
x3 = x3/x3(3);
x4 = x4/x4(3);

square = Rectangle2Square([x1 x2 x3 x4],W,L);     % ���ָ�����������ƽ���ϵ�����������������x1Ϊ����ԭ�㣩����ͶӰ����һ��ƽ����,�õ��������������

y1 = square(:,1);
y2 = square(:,2);
y3 = square(:,3);
y4 = square(:,4);

H2 = cross(V_vert,V_hori);   % ��������ƽ��ķ��߷��򣬴�ֱ��dx��dy�ķ���Ҳ����������ƽ�����ʧ��

V_diag1 = cross(cross(y1,y3),H2);   % �Խ����ϵ���ʧ��
V_diag2 = cross(cross(y2,y4),H2);

V_diag1 = V_diag1 / norm(V_diag1);  % ������ָ��Խ�����ʧ��ĵ�λʸ��
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

        % ֱ�ߺ͵����㷽��x'*l=0����η��̵Ľ�Ϊ����x'��С����ֵ��Ӧ��������ʸ����������x��������ʸ��
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


x1 = cross(hori_first,vert_first);   % ƽ��py_1��ƽ��px_1�Ľ��߷��򣬴�����ָ���һ���ǵ�
x2 = cross(hori_first,vert_last);    % ƽ��py_1��ƽ��px_N_x�Ľ��߷��򣬴�����ָ��x�������һ���ǵ�
x3 = cross(hori_last,vert_last);     % ƽ��py_N_y��ƽ��px_N_x�Ľ��߷��򣬴�����ָ��ԽǴ���Զ�Ľǵ�
x4 = cross(hori_last,vert_first);    % ƽ��py_N_y��ƽ��px_1�Ľ��߷��򣬴�����ָ��y�������һ���ǵ�

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
H2 = H2/norm(H2);        % ��ʧ�ߵĵ�λ���ꡣ������ĺ���ʧ�ߵ�ƽ��ΪPov����H2ΪPov�ĵ�λ��ʸ��

idx = abs(dx);
idx = find(idx==max(idx));

if V_hori(idx) < 0, V_hori = -V_hori; end;        % the point at the infinity of x axis
V_hori = sign(dx(idx))*V_hori;
V_hori = V_hori - dot(V_hori,H2)*H2;
V_hori = V_hori/norm(V_hori);
Rc = [V_hori, cross(H2,V_hori), H2];
omc = rodrigues(Rc);

% Deduce the translation vector: x1 is the origin
X1 = x1 / dot(H2,x1);    % d1 = dot(H2,x1)Ϊx1��ƽ��Pov�ľ��룬 X1 = x1/d1ʹ��X1��Pov�ľ���Ϊ1;
X2 = x2 / dot(H2,x2);    % d2 = dot(H2,x2)Ϊx2��ƽ��Pov�ľ��룬 X2 = x2/d2ʹ��X2��Pov�ľ���Ϊ1;
X3 = x3 / dot(H2,x3);    % d3 = dot(H2,x3)Ϊx3��ƽ��Pov�ľ��룬 X3 = x3/d3ʹ��X3��Pov�ľ���Ϊ1;
X4 = x4 / dot(H2,x4);    % d4 = dot(H2,x4)Ϊx4��ƽ��Pov�ľ��룬 X4 = x4/d4ʹ��X4��Pov�ľ���Ϊ1;

% ��Ȼx1,x2,x3,x4��X1,X2,X3,X4͸�Ӷ�Ӧ����X1,X2,X3,X4����ƽ��ƽ����Povƽ�棨ƽ��������ƽ�棩���Ҿ���Ϊ1��

scale = X1(3);     %   scale = x13/d1

X1 = X1/scale;      % �������ź�X1(3)=1,��X1�ڹ�һ��ƽ���ϣ���x1,x2,x3,x4���ڹ�һ��ƽ���ϣ���X1��x1�غϡ�
X2 = X2/scale;      % X1,X2,X3,X4����ƽ���X1��ƽ����Povƽ�棬���X1,X2,X3,X4��������������
X3 = X3/scale;      % һ���ı���x1,x2,x3,x4�ָ�Ϊƽ���ı���X1,X2,X3,X4
X4 = X4/scale;

u_hori = X2 - X1;  % ��������Ϊ����ʸ��
u_vert = X4 - X1;

scale = mean([W/norm(u_hori); L/norm(u_vert)]);
Tc = scale * X1;

return;






% N_x = n_sq_x+1;
% N_y = n_sq_y+1;
% W = abs(Xgrid(1,N_x)-Xgrid(1,1));               % X���������ܳ���ΪW, Y���������ܳ���ΪL
% L = abs(Xgrid(2,N_x*(N_y-1)+1)-Xgrid(2,1));
%
% [ideal_grid,V_hori,V_vert,H]=denoise_grid(grid_pts,n_sq_x,n_sq_y);
% x1 = [ideal_grid(:,1);1];                 % ƽ��py_1��ƽ��px_1�Ľ��߷��򣬴�����ָ���һ���ǵ�
% x2 = [ideal_grid(:,N_x);1];               % ƽ��py_1��ƽ��px_N_x�Ľ��߷��򣬴�����ָ��x�������һ���ǵ�
% x3 = [ideal_grid(:,end);1];               % ƽ��py_N_y��ƽ��px_N_x�Ľ��߷��򣬴�����ָ��ԽǴ���Զ�Ľǵ�
% x4 = [ideal_grid(:,N_x*(N_y-1)+1);1];     % ƽ��py_N_y��ƽ��px_1�Ľ��߷��򣬴�����ָ��y�������һ���ǵ�
