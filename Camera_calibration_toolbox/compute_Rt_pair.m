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
%          xpair: corresponding image points of two cameras: 2*npts*2 or 4*npts
%          fmat: Camera focal length of every camera: 2*2
%          cmat: Principal point coordinates of every Camera: 2*2
%          kmat: Di however,stortion coefficients of every camera: 5*2
%          alphavec: Skew coefficient of every camera: 1*2
%          hand: Handness of camera 2 wrt camera 1;
%
% OUTPUT:
%          Omc: the rotation vector (axis angle) from camera 1 to camera 2;
%          Tc: the translation vector from camera 1 to camera 2;
%          Emat: the essential matrix, x2n'*Emat*x1n=0;

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
npts = 20;
checkoutpic = 1;
maxfov = 90;        % unit: degree

X = [500*randn(2,npts); 2000+200*randn(1,npts)];      % unit: mm
Xm = mean(X,2); % aim of cameras
om = [zeros(3,1),randn(3,1)];
om(:,2) = om(:,2)/norm(om(:,2))*rand*pi;
hand = [1,sign(randn)]
Rt = rodrigues(om(:,2))';
if hand(2)~=1,
    Rt(3,:) = -Rt(3,:);
end;
direc = -Rt(:,3);    % the direction of -Z axis of camera 2
T = [zeros(3,1),-Rt'*(Xm+(2000+200*randn)*direc)];

% [x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha,dxdX] = project_points_mirror2(X,om,T,hand,f,c,k,alpha);
xxx = NaN(4,npts);
imageXY = 500+randi(500,2,2);
fov_angle = 25+randi(maxfov,1,2);
f = (imageXY/2)./repmat(tan(pi*fov_angle/360),2,1);
c = (imageXY-1)/2;
div = 3.^(repmat((1:5)',1,2))+5;
k = randn(5,2)./div;
alpha = 0.01*randn(1,2);
for i = 1:2,
    x = project_points_mirror2(X,om(:,i),T(:,i),hand(i),f(:,i),c(:,i),k(:,i),alpha(i));
    if checkoutpic,
        nx = imageXY(1,i);
        ny = imageXY(2,i);
        outind = x(1,:)>nx-0.5 | x(1,:)<-0.5 | x(2,:)>ny-0.5 | x(2,:)<-0.5;
        x(:,outind) = NaN;
    end;
    ii = (i-1)*2;
    xxx(ii+1:ii+2,:) = x;
end;

id = all(~isnan(xxx),1);
[om2, T2, E] = compute_Rt_pair(xxx(:,id),f,c,k,alpha,hand(2));
eom = om2-om(:,2)
s = T(:,2)./T2
