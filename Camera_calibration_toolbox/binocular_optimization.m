function [X,om2,T2,fc2,cc2,kc2,alpha2,estd,estd0] = binocular_optimization(xpair,om,T,hand,fc,cc,kc,alpha,...
                                                                           est_fc,center_optim,est_dist,est_alpha,est_aspect)
% BINOCULAR_OPTIMIZATION computes the optimized 3D structure X in a given binocular
% system, the parameters of the two cameras will be refined at the same time.
%
% INPUT:
%       xpair: corresponding image points of two cameras (2*npts*2 or 4*npts);
%       om: the initial rotation vector (axis angle) from camera 1 to camera 2 (3*1);
%       T: the initial translation vector from camera 1 to camera 2 (3*1);
%       hand: the handedness of camera 2 wrt camera 1 (1 or -1);
%       fc: the initial camera focal length of the two cameras (2*2);
%       cc: the initial principal point coordinates of the two cameras (2*2);
%       kc: the initial distortion coefficients of every camera (5*2);
%       alpha: the initial skew coefficient of every camera (1*2);
%       est_fc: switch to turn on/off the estimation of focal length (2*2);
%       center_optim: switch to turn on/off the estimation of principal point (1*2);
%       est_dist: switch to turn on/off the estimation of distortion coefficients (5*2);
%       est_alpha: switch to turn on/off the estimation of pixel skew (1*2);
%       est_aspect: switch to turn on/off the estimation of aspect ratio of focal length (1*2);
%
% OUTPUT:
%       X: the reconstructed 3d points in the 1st camera frame (3*npts)
%       om2: the refined rotation vector (axis angle) from camera 1 frame (3*1);
%       T2: the refined translation vector from camera 1 frame (3*1);
%       fc2: the refined camera focal length of the two cameras (2*2);
%       cc2: the refined principal point coordinates of the two cameras (2*2);
%       kc2: the refined distortion coefficients of every camera (5*2);
%       alpha2: the refined skew coefficient of every camera (1*2);
%
% Important functions called within that program:
% normalize_pixel: Computes the normalize image point coordinates.
%
% See also binocular_optimization2, compute_Rt_pair, compute_structure2, stereo_triangulation2.

% By ZPF @ZVR, 2017-8-24

MaxIter = 10; % Maximum number of iterations
bigeps = 1e-5;

if nargin<13,
    est_aspect = true(1,2);
    if nargin<12,
        est_alpha = true(1,2);
        if nargin<11,
            est_dist = true(5,2);
            if nargin<10,
                center_optim = true(1,2);
                if nargin<9,
                    est_fc = true(2,2);
                    if nargin<8,
                        alpha = [0,0];
                        if nargin<7,
                            kc = zeros(5,2);
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

[m,npts,n] = size(xpair);
if n==1,
    assert(m==4,'The 1st argument must contain corresponding points of two cameras!');
    X = xpair;
    xpair = zeros(2,npts,2);
    xpair(:,:,1) = X(1:2,:);
    xpair(:,:,2) = X(3:4,:);
elseif n==2,
    assert(m==2,'Unexpected dimension of the 1st argument!');
else
    disp('Unexpected dimension of the 1st argument!');
    return;
end;

assert(isequal(size(om),[3,1]), 'Unexpected dimension of the rotation vector!');
assert(isequal(size(T),[3,1]), 'Unexpected dimension of the translation vector!');
assert(isreal(hand) && isequal(abs(hand),1), 'Unexpected value of the handedness!');
assert(isequal(size(fc),[2,2]), 'Unexpected dimension of the 2 cameras'' focal length!');
assert(isequal(size(cc),[2,2]), 'Unexpected dimension of the 2 cameras'' principle points!');
assert(isequal(size(kc),[5,2]), 'Unexpected dimension of the 2 cameras'' distortion coefficients!');
assert(isequal(size(alpha),[1,2]), 'Unexpected dimension of the 2 cameras'' skew parameters!');
assert(isequal(size(est_fc),[2,2]), 'Unexpected dimension of the 9th argument!');
assert(isequal(size(center_optim),[1,2]), 'Unexpected dimension of the 10th argument!');
assert(isequal(size(est_dist),[5,2]), 'Unexpected dimension of the 11th argument!');
assert(isequal(size(est_alpha),[1,2]), 'Unexpected dimension of the 12th argument!');
assert(isequal(size(est_aspect),[1,2]), 'Unexpected dimension of the 13th argument!');

omc = [zeros(3,1),om];
T = T/norm(T);
Tc = [zeros(3,1),T];
handcc = [1,hand];
X = compute_structure2(xpair,omc,Tc,handcc,fc,cc,kc,alpha);
idx = all(~isnan(X),1);

param_up = reshape([fc; cc; alpha; kc; omc; Tc],32,1);
param = param_up;

% The following vector helps to select the variables to update:
ind = [est_fc; ones(2,1)*center_optim; est_alpha; est_dist; zeros(6,1), ones(6,1)];
ind(2,:) = ind(2,:).*(est_aspect | ~est_fc(1,:));
ind = logical(ind(:)');

% initial error before bundle adjustment
ex = []; % Global error vector
for pp = 1:2,
    x_kk = xpair(:,idx,pp);
    x = project_points_mirror2(X(:,idx),omc(:,pp),Tc(:,pp),handcc(pp),fc(:,pp),cc(:,pp),kc(:,pp),alpha(pp));
    ex_kk = x_kk - x;
    ex = [ex, ex_kk];
end;
if nargout > 8,
    estd0 = std(ex,0,2);
end;

lamda = 0.001; % set an initial value of the damping factor for the LM
updateJ = 1;
ex = ex(:);
ex2 = dot(ex,ex);
ex2 = ex2;
for iter = 1:MaxIter,
    % fprintf(1,'%d...',iter);
    if updateJ,
        % JJ2 = JJ'*JJ = U
        U = sparse([],[],[],32,32,512);
        ea = zeros(32,1);        % A'*ex

        % restore 3D points and camera parameters
        X = reshape(param,16,2);
        fc2 = X(1:2,:);
        cc2 = X(3:4,:);
        alpha2 = X(5,:);
        kc2 = X(6:10,:);
        omc = X(11:13,:);
        Tc = X(14:16,:);
        X = compute_structure2(xpair,omc,Tc,handcc,fc2,cc2,kc2,alpha2);
        for pp = 1:2,
            % load pixel points
            x_kk = xpair(:,idx,pp);
            if est_aspect(pp),
                [x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points_mirror2(X(:,idx),omc(:,pp),Tc(:,pp),handcc(pp),...
                                                                                fc2(:,pp),cc2(:,pp),kc2(:,pp),alpha2(pp));
            else
                [x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points_mirror2(X(:,idx),omc(:,pp),Tc(:,pp),handcc(pp),...
                                                                                fc2(1,pp),cc2(:,pp),kc2(:,pp),alpha2(pp));
                dxdf = repmat(dxdf,[1 2]);
            end;
            ex_kk = x_kk - x;
            Akk = [dxdf, dxdc, dxdalpha, dxdk, dxdom, dxdT];
            ii = (pp-1)*16;
            U(ii+1 : ii+16, ii+1 : ii+16) = Akk'*Akk;
            ea(ii+1 : ii+16) = Akk'*ex_kk(:);
        end;
        U = U(ind,ind);
        ea = ea(ind);
    end;
    U_lm = U + diag(lamda*diag(U));  % U + lamda*speye(size(U));
    param_innov = U_lm\ea;
    param_up(ind) = param(ind) + param_innov;     % updated parameters

    % compute reprojection error vector
    X = reshape(param_up,16,2);
    fc2 = X(1:2,:);
    cc2 = X(3:4,:);
    alpha2 = X(5,:);
    kc2 = X(6:10,:);
    omc = X(11:13,:);
    Tc = X(14:16,:);
    X = compute_structure2(xpair,omc,Tc,handcc,fc2,cc2,kc2,alpha2);
    ex = [];
    for pp = 1:2,
        % New intrinsic parameters:
        if ~est_aspect(pp) && all(est_fc(:,pp)),
            fc2(2,pp) = fc2(1,pp);
            ii = (pp-1)*16;
            intr_up(ii+2) = intr_up(ii+1);
        end;
        % load pixel points
        x_kk = xpair(:,idx,pp);
        x = project_points_mirror2(X(:,idx),omc(:,pp),Tc(:,pp),handcc(pp),fc2(:,pp),cc2(:,pp),kc2(:,pp),alpha2(pp));
        ex_kk = x_kk - x;
        ex = [ex, ex_kk];
    end;
    ex = ex(:);
    ex2_up = dot(ex,ex);
    if ex2_up < ex2,
        lamda = lamda/10;
        param = param_up;
        ex2 = ex2_up;
        updateJ = 1;
    else
        lamda = lamda*10;
        updateJ = 0;
    end;
    if norm(ea) < bigeps,
        break;
    end;
end;

%%%--------------------------- Computation of the error of estimation:

% fprintf(1,'\nEstimation of uncertainties...\n');
param = reshape(param, 16, 2);
fc2 = param(1:2,:);
cc2 = param(3:4,:);
alpha2 = param(5,:);
kc2 = param(6:10,:);
om2 = param(11:13,2);
T2 = param(14:16,2);
T2 = T2/norm(T2);
Tc = [zeros(3,1),T2];
omc = [zeros(3,1),om2];
X = compute_structure2(xpair,omc,Tc,handcc,fc2,cc2,kc2,alpha2,10);

% compute the reprojected errors
ex = [];
for pp = 2,
    x_kk = xpair(:,idx,pp);
    x = project_points_mirror2(X(:,idx),omc(:,pp),Tc(:,pp),handcc(pp),fc2(:,pp),cc2(:,pp),kc2(:,pp),alpha2(pp));
    ex_kk = x_kk - x;
    ex = [ex, ex_kk];
end;

if nargout > 7,
    estd = std(ex,0,2);
end;

return;

sigma_x = std(ex(:));

% Compute the the standard deviation of parameters from std(ex):
% ex = X-f(P),  cov(param,param) = inv((JJ'* inv(cov(ex,ex))* JJ))
JJ2_inv = eye(size(U))/U;

% Extraction of the final intrinsic paramaters:
param_error = zeros(32,1);           % take value as perfect if not optimized (no error).
param_error(ind) = 3*sqrt(abs(full(diag(JJ2_inv))))*sigma_x;     % 3 sigma principle
param_error = reshape(param_error,16,2);
ind = ~est_aspect & all(est_fc,1);
param_error(2,ind) = param_error(1,ind);

fc_error = param_error(1:2,:);
cc_error = param_error(3:4,:);
alpha_error = param_error(5,:);
kc_error = param_error(6:10,:);
om_error = param_error(11:13,2);
T_error = param_error(14:16,2);


return;



%% Test
np = 1000;
checkoutpic = 1;
flag = 0;
fov_angle = 70+randn*5;
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
imageXY = 500+randi(500,2,2);
f = [1;1].*((imageXY(1,:)/2)./tanfov_2)+10*randn(2,2);
c = (imageXY-1)/2+50*randn(2,2);
% k = zeros(5,2);
k = [randn(1,2)/20; randn(1,2)/100; randn(2,2)/200; randn(1,2)/1000];
alp = randi([-1,1],1,2).*rand(1,2)/10;
for i = 1:2,
    x = project_points_mirror2(X,om(:,i),T(:,i),hd(i),f(:,i),c(:,i),k(:,i),alp(i));
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
% [om2, T2] = compute_Rt_pair(xx,f,c,k,alp,hd(2));
f1 = f+100*randn;
c1 = (imageXY-1)/2;
k1 = zeros(5,2);
alp1 = zeros(1,2);
[om1, T1] = compute_Rt_pair(xx,f1,c1,k1,alp1,hd(2));
[XX,om2,T2,f2,c2,k2,alp2,e,e0] = binocular_optimization(xx,om1,T1,hd(2),f1,c1,k1,alp1,ones(2),ones(1,2),zeros(5,2),ones(1,2),ones(1,2));
XXlen = sqrt(sum(diff(XX,1,2).^2,1));
Xlen = sqrt(sum(diff(X,1,2).^2,1));
id = Xlen>1e-3;
s = mean(Xlen(id)./XXlen(id));
err = norm([reshape([f2-f;c2-c;alp2-alp;k2-k],20,1); om2-om(:,2); T2-T(:,2)/s]);
if err>1e-5,
    flag = 1;
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
        alpc = alp(pp);
        nx = imageXY(1,pp);
        ny = imageXY(2,pp);
        KK = [fc(1) fc(1)*alpc cc(1);0 fc(2) cc(2); 0 0 1];
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
