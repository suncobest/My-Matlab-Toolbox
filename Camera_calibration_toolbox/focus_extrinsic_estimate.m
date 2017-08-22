function [fc,omc,Tc] = focus_extrinsic_estimate(grid_pts,Xgrid,cc,f_ini,k_dist,two_focal,N_iter)

% Computes the focal length and the extrinsic parameters knowing the principal point and distortion factor k_dist
%
% grid_pts are the grid pixel coordinates. Xgrid are the world grid points.
% cc is the principal point. Set it to the center of image if you dont't know it.
% fc is the focal length being estimated. omc is the rotation of world plane, Tc is the translation.
%
% First estimate focus without distortion,then run iteration to compensate for distortion.
% The results will be less accurate when distortion is large.
%
% See also focus_extrinsic_estimate2, center_extrinsic_estimate, Distor2Calib.

% California Institute of Technology
% (c) Jean-Yves Bouguet - October 7th, 1997
% Beijing Institute of Technology
% Zhang Pengfei - May 4th, 2015


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

Xgrid = Xgrid(1:2,:);

flag = 0;
if norm(k_dist)~=0,
    flag = 1;
end;


% 一幅网格图像只能提供关于IAC的两个约束条件（一个平面上有两个虚圆点），因此只能估计焦距或主点，这里将图像中心作为主点估计焦距
fc = f_ini;

H = compute_homography_lm(grid_pts,Xgrid);

% matrix that subtract the principal point:
Sub_cc = [1 0 -cc(1);0 1 -cc(2);0 0 1];
H = diag([1/fc(1),1/fc(2),1]) * Sub_cc * H;

% 将图像坐标原点平移到主点，不考虑倾斜alpha_c，则有H=[h1,h2,h3] = diag(fc1,fc2,1)*[r1,r2,t]
% 此时内参数K=diag(fc1,fc2,1)；绝对二次曲线的像（IAC）为omega=inv(K')*inv(K)
% 齐次共轭点对u和v满足：u'*omega*v=0
% 4个点全部位于消失线上, h1'*diag([1/fc1^2,1/fc2^2,1])*h2 =0;
% h1'*diag([1/fc1^2,1/fc2^2,1])*h1 - h2'*diag([1/fc1^2,1/fc2^2,1])*h2 = (h1+h2)'*diag([1/fc1^2,1/fc2^2,1])*(h1-h2) = 0
% so h1+h2 and h1-h2 are conjugate points

invZc = mean([norm(H(:,1));norm(H(:,2))]);
H = H/invZc;

V_hori_pix = H(:,1);                  % h1
V_vert_pix = H(:,2);                  % h2
V_diag1_pix = H(:,1)+H(:,2);          % h1+h2
V_diag2_pix = H(:,1)-H(:,2);          % h1-h2

a1 = V_hori_pix(1);
b1 = V_hori_pix(2);
c1 = V_hori_pix(3);

a2 = V_vert_pix(1);
b2 = V_vert_pix(2);
c2 = V_vert_pix(3);

a3 = V_diag1_pix(1);
b3 = V_diag1_pix(2);
c3 = V_diag1_pix(3);

a4 = V_diag2_pix(1);
b4 = V_diag2_pix(2);
c4 = V_diag2_pix(3);

A = [a1*a2  b1*b2; a3*a4  b3*b4];
b = -[c1*c2;c3*c4];

if two_focal,
    f = sqrt(abs(1./(A\b)));
else
    f = sqrt(abs(1./(sum(A,2)\b))) * ones(2,1);
end;

H = diag([1/f(1),1/f(2),1])*H;
fc = fc.*f;

if flag,
    grid_centered = [grid_pts(1,:) - cc(1);grid_pts(2,:) - cc(2)];
    for i=1:N_iter,
        % Subtract principal point, and divide by the focal length:
        x_distort = [grid_centered(1,:)/fc(1); grid_centered(2,:)/fc(2)];
        x_undist = comp_distortion(x_distort,k_dist);
        H = compute_homography_lm(x_undist,Xgrid);                              % 多次迭代后H=s*[r1,r2,T]

        invZc = mean([norm(H(:,1));norm(H(:,2))]);
        H = H/invZc;

        V_hori_pix = H(:,1);                  % h1
        V_vert_pix = H(:,2);                  % h2
        V_diag1_pix = H(:,1)+H(:,2);          % h1+h2
        V_diag2_pix = H(:,1)-H(:,2);          % h1-h2

        a1 = V_hori_pix(1);
        b1 = V_hori_pix(2);
        c1 = V_hori_pix(3);

        a2 = V_vert_pix(1);
        b2 = V_vert_pix(2);
        c2 = V_vert_pix(3);

        a3 = V_diag1_pix(1);
        b3 = V_diag1_pix(2);
        c3 = V_diag1_pix(3);

        a4 = V_diag2_pix(1);
        b4 = V_diag2_pix(2);
        c4 = V_diag2_pix(3);

        A = [a1*a2  b1*b2; a3*a4  b3*b4];
        b = -[c1*c2;c3*c4];
        if two_focal,
            f = sqrt(abs(1./(A\b)));                            % 多次迭代后f=[1;1]
        else
            f = sqrt(abs(1./(sum(A,2)\b))) * ones(2,1);
        end;
        fc = fc.*f;
    end;
end;

invZc = mean([norm(H(:,1));norm(H(:,2))]);  % 求出h1，h2的平均长度，深度的倒数
H = H/invZc;
if H(9)<0,                                  % 世界坐标系原点的深度为正
    H=-H;
end;
h1 = H(:,1);
h1 = h1/norm(h1);
h2 = H(:,2)-dot(H(:,2),h1)*h1;
h2 = h2/norm(h2);
R = [h1 h2 cross(h1,h2)];
omc = rodrigues(R);
Tc = H(:,3);

return;



%% test
Nx = 5;
Ny = 5;
scale = randi(100);
[Ygrid,Xgrid]=meshgrid(((1:Ny)-1)*scale,((1:Nx)-1)*scale);
X0=[Xgrid(:)';Ygrid(:)';ones(1,Nx*Ny)];
K=[scale+10*rand, 0, scale/2+5*rand; 0, abs(scale-10*rand), scale/2-5*rand; 0 0 1]
T=[(randi(max(Nx,Ny),2,1)-(max(Nx,Ny)-1)/2);10]*scale;
R=rodrigues(randn(3,1)*1.5);
H=K*[R(:,1:2),T];
x=H*X0;
x=x./(ones(3,1)*x(3,:));
x1234=x(:,[1,Nx,Nx*Ny,Nx*(Ny-1)+1]);
[du,dv]=UnWarpPlane(x1234);
dot(du,dv)

kdist = [-0.1;0.02;-0.001;-0.002;-0.002];
y=add_distortion(x(1:2,:),kdist,[K(7);K(8)],[K(1);K(5)]);  % make sure center and fc is right.
y=y(1:2,:)+rand(2,Nx*Ny)/50;

xgrid = zeros(Nx,Ny);
ygrid = zeros(Nx,Ny);
xgrid(:) = x(1,:);
ygrid(:) = x(2,:);
figure(1);
plot(x(1,:),x(2,:),'r.')
hold on,axis equal
plot(xgrid,ygrid,'r-',xgrid',ygrid','r-')


xgrid(:) = y(1,:);
ygrid(:) = y(2,:);
figure(1);
plot(y(1,:),y(2,:),'g.')
hold on,axis equal
plot(xgrid,ygrid,'g-',xgrid',ygrid','g-')

set(gcf,'color','w');
axis tight

[fc,omc,Tc] = focus_extrinsic_estimate(y,X0,[K(7);K(8)],1,0,1);                      % do not consider distortion
[fc1,omc1,Tc1] = focus_extrinsic_estimate(y,X0,[K(7);K(8)],1,kdist,1);           % iterate giving the distortion
[fc2,omc2,Tc2] = focus_extrinsic_estimate2(y,X0,[K(7);K(8)],1,kdist,1);

err = [fc-[K(1);K(5)];reshape(rodrigues(omc),[],1)-R(:);Tc-T];
err = err'*err
err1 = [fc1-[K(1);K(5)];reshape(rodrigues(omc1),[],1)-R(:);Tc1-T];
err1 = err1'*err1
err2 = [fc2-[K(1);K(5)];reshape(rodrigues(omc2),[],1)-R(:);Tc2-T];
err2 = err2'*err2
