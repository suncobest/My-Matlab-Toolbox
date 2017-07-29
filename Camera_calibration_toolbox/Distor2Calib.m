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
x_grid(:) = grid_pts_centered(1,:);   % ����ƽ������Ϊԭ�������ϵ
y_grid(:) = grid_pts_centered(2,:);

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
    [~,~,V] = svd(X);    % ��ƽ������ȷ���ƽ��1�����࿴��1��ԭ��ƽ������Ϊԭ�㣬��Ϊ���ģ���ʼλ�ã�
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

% ֱ�ߺ͵����㷽��x'*l=0����η��̵Ľ�Ϊ����x'��С����ֵ��Ӧ��������ʸ����������x��������ʸ��
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


