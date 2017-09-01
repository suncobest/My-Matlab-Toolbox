function [X,om2,T2,estd,estd0] = binocular_optim1D_extrinsic(xpair,rodlen,om,T,hand,fc,cc,kc,alpha)
% BINOCULAR_OPTIM1D_EXTRINSIC computes the optimized 3D structure X in a given binocular
% system, the parameters of the two cameras will be refined at the same time.
%
% INPUT:
%       xpair: corresponding image points of two cameras (2*npts*2 or 4*npts);
%       np1D: number of points on the rod;
%       rodlen: the length of the rod;
%       om: the initial rotation vector (axis angle) from camera 1 to camera 2 (3*1);
%       T: the initial translation vector from camera 1 to camera 2 (3*1);
%       hand: the handedness of camera 2 wrt camera 1 (1 or -1);
%       fc: the initial camera focal length of the two cameras (2*2);
%       cc: the initial principal point coordinates of the two cameras (2*2);
%       kc: the initial distortion coefficients of every camera (5*2);
%       alpha: the initial skew coefficient of every camera (1*2);
%
% OUTPUT:
%       X: the reconstructed 3d points in the 1st camera frame (3*npts);
%       om2: the refined rotation vector (axis angle) from camera 1 frame (3*1);
%       T2: the refined translation vector from camera 1 frame (3*1);
%
% Important functions called within that program:
% normalize_pixel: Computes the normalize image point coordinates.
%
% See also binocular_1D_optim, binocular_optimization2, compute_Rt_pair, compute_structure2, stereo_triangulation2.

% By ZPF @ZVR, 2017-8-24

MaxIter = 10; % Maximum number of iterations
bigeps = 1e-5;
thresh_cond = 1e10;

if nargin<9,
    alpha = [0,0];
    if nargin<8,
        kc = zeros(5,2);
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

np1D = length(rodlen);
assert(np1D>1 && ~mod(npts,np1D),'Unexpected number of points on the rod!');
assert(isequal(size(om),[3,1]), 'Unexpected dimension of the rotation vector!');
assert(isequal(size(T),[3,1]), 'Unexpected dimension of the translation vector!');
assert(isreal(hand) && isequal(abs(hand),1), 'Unexpected value of the handedness!');
assert(isequal(size(fc),[2,2]), 'Unexpected dimension of the 2 cameras'' focal length!');
assert(isequal(size(cc),[2,2]), 'Unexpected dimension of the 2 cameras'' principle points!');
assert(isequal(size(kc),[5,2]), 'Unexpected dimension of the 2 cameras'' distortion coefficients!');
assert(isequal(size(alpha),[1,2]), 'Unexpected dimension of the 2 cameras'' skew intreters!');

omc = [zeros(3,1),om];
T = T/norm(T);
Tc = [zeros(3,1),T];
handcc = [1,hand];
X = compute_structure2(xpair,omc,Tc,handcc,fc,cc,kc,alpha);

% solve the origin and direction of rods
nima = npts/np1D;
X = reshape(X,[3,np1D,nima]);
ind = reshape(all(all(~isnan(X),1),2),1,nima);
X = X(:,:,ind);
Xlen = permute(sqrt(sum(diff(X,[],2).^2,1)),[3,2,1]);  % length of rod (with ||T2||=1)
s = mean(diff(rodlen)./mean(Xlen,1));
T = T*s;
X = X*s;
Tc(:,2) = T;
Xo = X(:,1:np1D:end);
Xn = X(:,np1D:np1D:end)-Xo;
Xn = Xn./(ones(3,1)*sqrt(sum(Xn.^2,1)));
thph = cartesian2spherical(Xn);
thph = thph(2:3,:);
X = reshape(permute(Xo+Xn.*reshape(rodlen,[1,1,np1D]), [1,3,2]), 3,npts);

idx = reshape(ind(ones(np1D,1),:),1,npts);
nact = sum(ind);
nact5 = 5*nact;
nid2 = 2*nact*np1D;

intr_up = [om; T];
extr_up = reshape([Xo; thph],nact5,1);
intr = intr_up;
extr = extr_up;

% initial error before bundle adjustment
ex = []; % Global error vector
for pp = 1:2,
    x_kk = xpair(:,idx,pp);
    x = project_points_mirror2(X,omc(:,pp),Tc(:,pp),handcc(pp),fc(:,pp),cc(:,pp),kc(:,pp),alpha(pp));
    ex_kk = x_kk - x;
    ex = [ex, ex_kk];
end;
if nargout > 4,
    estd0 = std(ex,0,2);
end;

lamda = 0.001; % set an initial value of the damping factor for the LM
updateJ = 1;
ex = ex(:);
ex2 = dot(ex,ex);
for iter = 1:MaxIter,
    % fprintf(1,'%d...',iter);
    if updateJ,
        % JJ2 = JJ'*JJ = [U, W; W', V]
        U = zeros(6);
        V = sparse([],[],[],nact5,nact5,5*nact5);
        W = sparse([],[],[],6,nact5,6*nact5);
        ea = zeros(6,1);        % A'*ex
        eb = zeros(nact5,1);    % B'*ex

        % restore 3D points and camera parameters
        omc = [zeros(3,1),intr(1:3)];
        Tc = [zeros(3,1),intr(4:6)];
        X = reshape(extr,5,nact);
        Xo = X(1:3,:);
        thph = X(4:5,:);
        [X,dXdXo,dXdtp] = gen_1D_points(Xo,thph,rodlen);
        for pp = 1:2,
            % load pixel points
            x_kk = xpair(:,idx,pp);
            [x,dxdom,dxdT,~,~,~,~,dxdX] = project_points_mirror2(X,omc(:,pp),Tc(:,pp),handcc(pp),fc(:,pp),cc(:,pp),kc(:,pp),alpha(pp));
            ex_kk = x_kk - x;
            Bkk = sparse([],[],[],nid2,nact5);
            dxdXo = dxdX*dXdXo;
            dxdtp = dxdX*dXdtp;
            for jj = 1:nact,
                Bkk(:,(jj-1)*5+1:jj*5) = [dxdXo(:,(jj-1)*3+1:jj*3),dxdtp(:,(jj-1)*2+1:jj*2)];
            end;
            V = V + Bkk'*Bkk;
            eb = eb + Bkk'*ex_kk(:);
            if pp == 2,
                Akk = [dxdom, dxdT];
                U = Akk'*Akk;
                W = Akk'*Bkk;
                ea = Akk'*ex_kk(:);
            end;
        end;
    end;
    U_lm = U + diag(lamda*diag(U));  % U + lamda*speye(size(U));
    V_lm = V + diag(lamda*diag(V));  %  V + lamda*speye(size(V));
    Y = W/V_lm;

    intr_innov = (U_lm-Y*W')\(ea-Y*eb);              % da
    extr_innov = V_lm\(eb-W'*intr_innov);            % db
    intr_up = intr + intr_innov;     % updated parameters
    extr_up = extr + extr_innov;

    % compute reprojection error vector
    omc = [zeros(3,1),intr_up(1:3)];
    Tc = [zeros(3,1),intr_up(4:6)];
    X = reshape(extr_up,5,nact);
    X = gen_1D_points(X(1:3,:),X(4:5,:),rodlen);
    ex = [];
    for pp = 1:2,
        % load pixel points
        x_kk = xpair(:,idx,pp);
        x = project_points_mirror2(X,omc(:,pp),Tc(:,pp),handcc(pp),fc(:,pp),cc(:,pp),kc(:,pp),alpha(pp));
        ex_kk = x_kk - x;
        ex = [ex, ex_kk];
    end;
    ex = ex(:);
    ex2_up = dot(ex,ex);
    if ex2_up < ex2,
        lamda = lamda/10;
        intr = intr_up;
        extr = extr_up;
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
om2 = intr(1:3);
T2 = intr(4:6);
Tc = [zeros(3,1),T2];
omc = [zeros(3,1),om2];
if condest(V_lm)>thresh_cond,
    X = compute_structure2(xpair,omc,Tc,handcc,fc,cc,kc,alpha,10);
else
    X = reshape(extr,5,nact);
    Xo = X(1:3,:);
    thph = X(4:5,:);
    X = NaN(3,npts);
    X(:,idx) = gen_1D_points(Xo,thph,rodlen);
end;

% compute the reprojected errors
ex = [];
for pp = 2,
    x_kk = xpair(:,idx,pp);
    x = project_points_mirror2(X(:,idx),omc(:,pp),Tc(:,pp),handcc(pp),fc(:,pp),cc(:,pp),kc(:,pp),alpha(pp));
    ex_kk = x_kk - x;
    ex = [ex, ex_kk];
end;

if nargout > 3,
    estd = std(ex,0,2);
end;


return;



%% Test
nim = 1000;
np1d = 1+randi(4);
np = nim*np1d;
rlen = cumsum(0:np1d-1)*50;
l = 3000;
fov_angle = 70;
tanfov_2 = tan(pi*fov_angle/360);
Xori = l*tanfov_2/2*(2*rand(3,nim)-ones(3,1))+[0;0;l];  % the origin of rods
Xdir = randn(3,nim);   % the direction of rods
Xdir = Xdir./(ones(3,1)*sqrt(sum(Xdir.^2,1)));
Xrod = reshape(permute(Xori+Xdir.*reshape(rlen,[1,1,np1d]), [1,3,2]), 3,[]);

checkoutpic = 1;
flag = 0;
Xm = mean(Xrod,2); % aim of cameras
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
% k = zeros(5,2);
k = [randn(1,2)/20; randn(1,2)/100; randn(2,2)/200; randn(1,2)/1000];
alp = randi([-1,1],1,2).*rand(1,2)/10;
for i = 1:2,
    x = project_points_mirror2(Xrod,om(:,i),T(:,i),hd(i),f(:,i),c(:,i),k(:,i),alp(i));
    if checkoutpic,
        nx = imageXY(1,i);
        ny = imageXY(2,i);
        outid = x(1,:)>nx-0.5 | x(1,:)<-0.5 | x(2,:)>ny-0.5 | x(2,:)<-0.5;
        x(:,outid) = NaN;
    end;
    xx(1:2,:,i) = x;
end;

id = all(reshape(all(all(~isnan(xx),3),1),np1d,nim),1);
id = reshape(id(ones(np1d,1),:),1,np);
xx = xx(:,id,:);
Xrod = Xrod(:,id);
f1 = f+100*randn;
c1 = (imageXY-1)/2;
k1 = zeros(5,2);
alp1 = zeros(1,2);
[om1, T1] = compute_Rt_pair(xx,f1,c1,k1,alp1,hd(2));
[XX,om2,T2,e,e0] = binocular_optim1D_extrinsic(xx,rlen,om1,T1,hd(2),f1,c1);
XXlen = sqrt(sum(diff(XX,1,2).^2,1));
Xlen = sqrt(sum(diff(Xrod,1,2).^2,1));
id = Xlen>1e-3;
s = mean(Xlen(id)./XXlen(id));
err = norm([om2-om(:,2); T2-T(:,2)]);
if err>1e-1,
    flag = 1;
end;

if flag,
    delta = 0.2;        % shift ot text
    nu = l/10;
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
    plot3(Xrod(1,:),Xrod(3,:),-Xrod(2,:), '+');
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
