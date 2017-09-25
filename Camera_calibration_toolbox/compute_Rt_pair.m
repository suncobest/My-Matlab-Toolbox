function [Omc,Tc,Emat] = compute_Rt_pair(xpair,fmat,cmat,kmat,alphavec,hand)
% COMPUTE_RT_PAIR - computes essential matrix from corresponding points of two cameras
% and decompose it into R and t. Where
%         Xc2 = R*hand*Xc1+t.
% Function computes the essential matrix from 8 or more matching points in
% a stereo pair of images.  The normalised 8 point algorithm given by
% Hartley and Zisserman p265 is used.
% Decompose essential matrix into Rt: p257-p259.
% Linear triangulation methods: p312.
%
% Arguments:
%          xpair: corresponding image points of two cameras (2*npts*2 or 4*npts);
%          fmat: the camera focal length of every camera (2*2);
%          cmat: the principal point coordinates of every Camera (2*2);
%          kmat: the distortion coefficients of every camera (5*2);
%          alphavec: the skew coefficient of every camera (1*2);
%          hand: the handedness of camera 2 wrt camera 1 (1 or -1);
%
% OUTPUT:
%          Omc: the rotation vector (axis angle) from camera 1 to camera 2;
%          Tc: the translation vector from camera 1 to camera 2;
%          Emat: the essential matrix, x2n'*Emat*x1n=0;
%
% Method:
%         First normalizing all point's coordinates;
%         using the 8-point correspondence algorithm to compute the essential matrix;
%         decomposing the essential matrix into rotation and translation;
%         using DLT algorithm to compute 3D structure and choose the right solution.

% By ZPF @ZVR, 2017-8-4

[m,npts,n] = size(xpair);
if n==1,
    assert(m==4,'The 1st argument must contain corresponding points of two cameras!');
    x1 = xpair(1:2,:);
    x2 = xpair(3:4,:);
elseif n==2,
    assert(m==2,'Unexpected dimension of the 1st argument!');
    x1 = xpair(:,:,1);
    x2 = xpair(:,:,2);
else
    disp('Unexpected dimension of the 1st argument!');
    return;
end;

if npts < 8,
    error('At least 8 points are needed to compute the fundamental matrix!');
end;

[m, n] = size(fmat);
if m==2 && n==1,
    fmat = fmat(:,ones(1,2));
else
    assert(m==2 && n==2,'Unexpected dimension of the focal length matrix!');
end;

[m, n] = size(cmat);
if m==2 && n==1,
    cmat = cmat(:,ones(1,2));
else
    assert(m==2 && n==2,'Unexpected dimension of the principal point matrix!');
end;

if nargin < 6,
    hand = 1;
    if nargin < 5,
        alphavec = zeros(1,2);
        if nargin < 4,
            kmat = zeros(5,2);
        end;
    end;
end;

assert(abs(hand)==1,'The relative handedness is assumed to be 1 or -1!');

[m, n] = size(alphavec);
if m==1 && n==1,
    alphavec = alphavec(ones(1,2));
else
    assert(m==1 && n==2,'Unexpected dimension of aspect ratio vector!');
end;

[m, n] = size(kmat);
if m==5 && n==1,
    kmat = kmat(:,ones(1,2));
else
    assert(m==5 && n==2,'Unexpected dimension of the distortion matrix!');
end;

x1n = normalize_pixel(x1,fmat(:,1),cmat(:,1),kmat(:,1),alphavec(1));
x2n = normalize_pixel(x2,fmat(:,2),cmat(:,2),kmat(:,2),alphavec(2));

% Build the constraint matrix
A = [x2n(1,:)'.*x1n(1,:)',  x2n(1,:)'.*x1n(2,:)',  x2n(1,:)', ...
     x2n(2,:)'.*x1n(1,:)',  x2n(2,:)'.*x1n(2,:)',  x2n(2,:)', ...
     x1n(1,:)',             x1n(2,:)',           ones(npts,1)];

[~,~,V] = svd(A,0);
Emat = reshape(V(:,9),3,3)';

% Enforce constraint that fundamental matrix has rank 2 by performing
% a svd and then reconstructing with the two largest singular values.
[U,~,V] = svd(Emat);
Emat = U*diag([1, 1, 0])*V';

% decompose essential matrix into R and T: ||T|| = 1.
W = [0, -1, 0;
     1, 0, 0;
     0, 0, 1];
R = reshape([U*W*V', U*W'*V', U*W*V', U*W'*V'],[3,3,4]);
RT = [R,reshape(U(:,3)*[1, 1, -1, -1],[3,1,4])];
if hand~=1, % det(R)==-1
    R(:,3,:) = -R(:,3,:);
end;
if det(R(:,:,1))<0,
    R = -R;
    RT = -RT;
end;

x1 = x1n(:,1);
x2 = x2n(:,1);
id = 0;
% linear triangulation of two views
for pp=1:4,
    A = zeros(4);
    A(1:2,:) = x1*[0,0,1,0]-[1,0,0,0; 0,1,0,0];
    A(3:4,:) = x2*RT(3,:,pp)-RT(1:2,:,pp);
    [~,~,V] = svd(A);
    X = V(1:3,4)/V(4,4);
    Y = RT(:,:,pp)*[X;1];
    if all([X(3),Y(3)]>0),
        id = pp;
        break;
    end;
end;

Omc = rodrigues(R(:,:,id));
Tc = RT(:,4,id);

return;



%% Test
np = 500;
checkoutpic = 1;
div = [10; 50; 100; 100; 500]*[1,1];
flag = 0;
for count = 1:1000,
    fov_angle = 60+randn*10;
    tanfov_2 = tan(pi*fov_angle/360);
    l = 4000+randi([-1000,1000]);
    X = l*tanfov_2/2*(2*rand(3,np)-ones(3,1))+[0;0;l];      % unit: mm
    Xm = mean(X,2); % aim of cameras
    theta = pi/(rand*10+1);
    om = [zeros(3,1),theta*[cos(theta);0;-sin(theta)]];
    hd = [1,sign(randn)];
    Rt = rodrigues(om(:,2))';
    if hd(2)~=1,
        Rt(3,:) = -Rt(3,:);
    end;
    direc = -Rt(:,3);    % the direction of -Z axis of camera 2
    T = [zeros(3,1),-Rt'*(Xm+l*direc)];

    xx = NaN(2,np,2);
    imageXY = 800+randi(500,2,2);
    f = [1;1].*((imageXY(1,:)/2)./tanfov_2)+10*randn(2,2);
    c = (imageXY-1)/2+50*randn(2,2);
    k = randn(5,2)./div;
    % k = zeros(5,2);

    alpha = 0.01*randn(1,2);
    for i = 1:2,
        x = project_points_mirror2(X,om(:,i),T(:,i),hd(i),f(:,i),c(:,i),k(:,i),alpha(i));
        if checkoutpic,
            nx = imageXY(1,i);
            ny = imageXY(2,i);
            outind = x(1,:)>nx-0.5 | x(1,:)<-0.5 | x(2,:)>ny-0.5 | x(2,:)<-0.5;
            x(:,outind) = NaN;
        end;
        xx(1:2,:,i) = x;
    end;

    id = all(all(~isnan(xx),3),1);
    xx = xx(:,id,:);
    X = X(:,id);
    [om2, T2] = compute_Rt_pair(xx,f,c,k,alpha,hd(2));
    XX = compute_structure2(xx,[zeros(3,1),om2],[zeros(3,1),T2],hd,f,c,k,alpha);
    XXlen = sqrt(sum(diff(XX,1,2).^2,1));
    Xlen = sqrt(sum(diff(X,1,2).^2,1));
    id = Xlen>1e-3;
    s = mean(Xlen(id)./XXlen(id));
    err = norm([om2;T2]-[om(:,2);T(:,2)/s]);
    if err>1e-1,
        flag = 1;
        break;
    end;
end;

if flag,
    delta = 0.2;        % shift ot text
    nu = s/10;
    figure(4);
    clf;
    hold on;
    BASE = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]*nu;  % 原点 x基矢量 原点 y基矢量 原点 z基矢量
    DeltaXYZ = [BASE(:,[2,4,6])*(1+delta), -[1;1;1]*nu*delta];     % text位置: x, y, z, o
    for pp = 1:2,
        fc = f(:,pp);
        cc = c(:,pp);
        alpha_c = alpha(pp);
        nx = imageXY(1,pp);
        ny = imageXY(2,pp);
        KK = [fc(1) fc(1)*alpha_c cc(1);0 fc(2) cc(2); 0 0 1];
        IP = KK\[0 nx-1 nx-1 0 0; 0 0 ny-1 ny-1 0; 1 1 1 1 1]*nu;
        IP = reshape([IP;BASE(:,1)*ones(1,5);IP],3,15);
        % Change of reference: wrt reference frame
        Rckk = rodrigues(om(:,pp));
        if hd(pp)~=1,
            Rckk(:,3)=-Rckk(:,3);
        end;
        Twkk = T(:,pp);
        BASE_kk = Rckk'*(BASE-Twkk(:,ones(1,6)));       % Xw=Rk'*(Xk-Tk)
        DELTA_kk = Rckk'*(DeltaXYZ-Twkk(:,ones(1,4)));
        IP_kk = Rckk'*(IP -Twkk(:,ones(1,15)));
        figure(4);
        plot3(BASE_kk(1,:),BASE_kk(3,:),-BASE_kk(2,:),'b-','linewidth',1);
        plot3(IP_kk(1,:),IP_kk(3,:),-IP_kk(2,:),'r-','linewidth',1);
        text(DELTA_kk(1,1),DELTA_kk(3,1),-DELTA_kk(2,1),'X','HorizontalAlignment','center');
        text(DELTA_kk(1,2),DELTA_kk(3,2),-DELTA_kk(2,2),'Y','HorizontalAlignment','center');
        text(DELTA_kk(1,3),DELTA_kk(3,3),-DELTA_kk(2,3),'Z','HorizontalAlignment','center');
        text(DELTA_kk(1,4),DELTA_kk(3,4),-DELTA_kk(2,4),['Cam-' num2str(pp)],'HorizontalAlignment','center');
    end;
    plot3(X(1,:),X(3,:),-X(2,:), '+');
    az = 50;
    el = 20;
    figure(4);
    rotate3d on;
    grid on;
    title('Extrinsic results of the camera system');
    set(4,'color',[1 1 1],'Name','3D','NumberTitle','off');
    axis equal vis3d tight;
    xlabel('X');
    ylabel('Z');
    zlabel('-Y');
    view(az,el);
    hold off;
    err
end;
