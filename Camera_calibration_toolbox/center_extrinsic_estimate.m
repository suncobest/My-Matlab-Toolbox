function [cc,omc,Tc] = center_extrinsic_estimate(grid_pts,Xgrid,fc,c_ini,k_dist,N_iter)

% Computes the principal point and the extrinsic parameters giving the
% focal length and distortion factor k_dist
%
% grid_pts are the grid pixel coordinates. Xgrid are the world grid points.
% fc is the focal length being given. You can set the initial principal point c_ini
% to the center of image.
%
% cc is the principal point being estimated.
% omc is the rotation of world plane, Tc is the translation.
%
% Fisrt estimate center without distortion,then run iteration to compensate for distortion.
% The results will be less accurate when distortion is large. 
%
% See also focus_extrinsic_estimate, focus_extrinsic_estimate2.

% Beijing Institute of Technology
% Zhang Pengfei - May 9th, 2015

if nargin < 6,
    N_iter = 10;
    if nargin < 5,
        k_dist = 0;
        if nargin < 4,
            c_ini = [0;0];
        end;
    end;
end;

[mx,nx] = size(grid_pts);
[mX,nX] = size(Xgrid);
assert(mx==2||mx==3,'Unexpected points dimension!');
assert(mX==2||mX==3,'Unexpected points dimension!');
assert(nx==nX,'Number of 2D and 3D points unmatched!');
assert(length(c_ini)==2,'Unexpected principal point dimension!');

if mx == 3,
    grid_pts = grid_pts./(ones(3,1)*grid_pts(3,:));
end;

Xgrid = Xgrid(1:2,:);

if length(fc)==1,
    fc=[fc;fc];
end;

if length(k_dist)>1,
    k = zeros(5,1);
    k(k_dist~=0) = k_dist(k_dist~=0);
    k_dist = k;
    if norm(k(2:5))~=0,
        fprintf(1,'\nNote: High order of distortion coefficient detected! Make sure the distortion coefficient is right!\n');
    end;
end;

flag = 0;
if norm(k_dist)~=0,
    flag = 1;
end;


cc = c_ini;
H = compute_homography_lm(grid_pts,Xgrid);              % H = [h1,h2,h3]=sK(r1,r2,t)

% matrix that subtract the principal point:
Sub_cc = [1 0 -cc(1);0 1 -cc(2);0 0 1];
H = diag([1/fc(1),1/fc(2),1]) * Sub_cc * H;

% 一幅网格图像(或平行网格的图像)只能提供关于IAC的两个约束条件（一个平面上有两个虚圆点），因此只能估计焦距或主点，这里已知焦距估计主点
% h1'*IAC*h2=0; h1'*IAC*h1=h2'*IAC*h2; IAC=inv(K')*inv(K), inv(K)= [1/f1 0 -c1/f1;0 1/f2 -c2/f2;0 0 1];
% IAC=[w1,0,v1;0,w2,v2;v1,v2,v3], w1=1/f1^2, w2=1/f2^2, v1=-c1/f1^2, v2=-c2/f2^2, v3=(c1/f1)^2+(c2/f2)^2+1
% (h11*h32+h31*h12)*c1/f1^2+(h21*h32+h31*h22)*c2/f2^2 -h31*h32*((c1/f1)^2+(c2/f2)^2)= h11*h12/f1^2+h21*h22/f2^2+h31*h32,
% 2*(h11*h31-h12*h32)*c1/f1^2+2*(h21*h31-h22*h32)*c2/f2^2 -(h31^2-h32^2)*((c1/f1)^2+(c2/f2)^2)= (h11^2-h12^2)/f1^2+(h21^2-h22^2)/f2^2+h31^2-h32^2;
%
% 迭代时，K=[1 0 c1/f1;0 1 c2/f2;0 0 1], IAC=[1,0,v1;0,1,v2;v1,v2,v3],v1=-c1/f1, v2=-c2/f2, v3=(c1/f1)^2+(c2/f2)^2+1;
% 令x=c1/f1,y=c2/f2; h1'*IAC*h2=0; h1'*IAC*h1=h2'*IAC*h2;
% (h11*h32+h31*h12)*x+(h21*h32+h31*h22)*y = h1'*h2 + h31*h32*(x^2+y^2)
% 2*(h11*h31-h12*h32)*x+2*(h21*h31-h22*h32)*y = h1'*h1-h2'*h2 + (h31^2-h32^2)*((x^2+y^2);

h1 = H(:,1);
h2 = H(:,2);
A = [h1(1)*h2(3)+h1(3)*h2(1), h1(2)*h2(3)+h1(3)*h2(2), -h1(3)*h2(3); 2*(h1(1)*h1(3)-h2(1)*h2(3)), 2*(h1(2)*h1(3)-h2(2)*h2(3)), h2(3)^2-h1(3)^2];
b = [h1'*h2; h1'*h1-h2'*h2; 0];
c = [0;0;0];                                   % c=[x;y;x^2+y^2]
for j=1:10,
    AA = [A;c(1),c(2),-1];
    c = AA\b;
end;
c_over_f = c(1:2);

Sub_cc = [1 0 -c_over_f(1);0 1 -c_over_f(2);0 0 1];
H = Sub_cc * H;
cc = c_over_f.*fc + cc;

if flag,
    for i=1:N_iter,
        % Subtract principal point, and divide by the focal length:
        x_distort = [(grid_pts(1,:) - cc(1))/fc(1); (grid_pts(2,:) - cc(2))/fc(2)];
        x_undist = comp_distortion(x_distort,k_dist);
        H = compute_homography_lm(x_undist,Xgrid);                                      % 多次迭代后H=s*[r1,r2,T]
        
        h1 = H(:,1);
        h2 = H(:,2);
        A = [h1(1)*h2(3)+h1(3)*h2(1), h1(2)*h2(3)+h1(3)*h2(2), -h1(3)*h2(3); 2*(h1(1)*h1(3)-h2(1)*h2(3)), 2*(h1(2)*h1(3)-h2(2)*h2(3)), h2(3)^2-h1(3)^2];
        b = [h1'*h2; h1'*h1-h2'*h2; 0];
        c = [0;0;0];                                   % c=[x;y;x^2+y^2]
        for j=1:10,
            AA = [A;c(1),c(2),-1];
            c = AA\b;
        end;
        c_over_f = c(1:2);                             % 多次迭代后c_over_f=[0;0]       
        cc = c_over_f.*fc + cc;
    end;
end;

invZc = mean([norm(H(:,1));norm(H(:,2))]);  % 求出h1，h2的平均长度，深度的倒数
H = H/invZc;
if H(9)<0,
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
Nx = 6;
Ny = 5;
scale = randi(100);
[Ygrid,Xgrid]=meshgrid(((1:Ny)-1)*scale,((1:Nx)-1)*scale);
X0=[Xgrid(:)';Ygrid(:)';ones(1,Nx*Ny)];
K=[scale+10*rand, 0, scale/2+5*rand; 0, abs(scale-10*rand), scale/2-5*rand; 0 0 1];
T=[(randi(max(Nx,Ny),2,1)-(max(Nx,Ny)-1)/2);14]*scale;
R=rodrigues(randn(3,1)*0.3);
H=K*[R(:,1:2),T];
x=H*X0;
x=x./(ones(3,1)*x(3,:));
x1234=x(:,[1,Nx,Nx*Ny,Nx*(Ny-1)+1]);
[du,dv]=UnWarpPlane(x1234);
dot(du,dv)

kdist = [-0.15;0.01;0.02;-0.01;0.002];
y=add_distortion(x(1:2,:),kdist,[K(7);K(8)],[K(1);K(5)]);  % make sure center and fc is right.
y=y(1:2,:)+rand(2,Nx*Ny)/50;

xgrid = zeros(Nx,Ny);
ygrid = zeros(Nx,Ny);
xgrid(:) = x(1,:);
ygrid(:) = x(2,:);
figure(1);
plot(x(1,:),x(2,:),'r.')                               % accurate data
hold on,axis equal
plot(xgrid,ygrid,'r-',xgrid',ygrid','r-')


xgrid(:) = y(1,:);
ygrid(:) = y(2,:);
figure(1);
plot(y(1,:),y(2,:),'g.')                              % ditorted data
hold on,axis equal
plot(xgrid,ygrid,'g-',xgrid',ygrid','g-')

set(gcf,'color','w');
axis tight

[cc,omc,Tc] = center_extrinsic_estimate(y,X0,[K(1);K(5)]);
[cc1,omc1,Tc1] = center_extrinsic_estimate(y,X0,[K(1);K(5)],[0;0],kdist);
err = [cc-[K(7);K(8)];reshape(rodrigues(omc),[],1)-R(:);Tc-T];
err = err'*err
err1 = [cc1-[K(7);K(8)];reshape(rodrigues(omc1),[],1)-R(:);Tc1-T];
err1 = err1'*err1



