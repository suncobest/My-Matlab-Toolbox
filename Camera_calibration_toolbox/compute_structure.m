function [Xmat, goodpts] = compute_structure(xbody,Qmat,Tmat,handvec,fmat,cmat,kmat,alphavec,MaxIter)
% compute_structure
%
% Computes the 3D structure Xmat in a back projective way, given corresponding pixel
% points and parameters of n_cam cameras.
%
% INPUT: xbody: corresponding pixel feature stored in n_cam page of a cuboid: 2*npts*n_cam
%       Missing pixel data is denoted as NaN;
%       Qmat: the quaternion orientation of world frame in every camera frame: 4*n_cam
%       Tmat: the position of world frame origin in every camera frame: 3*n_cam
%       fmat: Camera focal length of every camera: 2*n_cam
%       cmat: Principal point coordinates of every Camera: 2*n_cam
%       kmat: Distortion coefficients of every camera: 5*n_cam
%       alphavec: Skew coefficient of every camera: 1*n_cam
%       handvec: Handness of every camera frame wrt world reference frame (1 or -1);
%
% OUTPUT:   Xmat: the restored 3d structure in world reference frame: 3*npts
%
% Method: First computes the normalized point coordinates, then computes the 3D
% structure using DLT algorithm.
%
% Important functions called within that program:
% normalize_pixel: Computes the normalize image point coordinates.
% trans_quat_mat: transform rotation matrix to quternion vector, or vice versa.
%
% See also compute_structure2, stereo_triangulation, stereo_triangulation2.

[two,npts,n_cam] = size(xbody);
assert(two==2,'Dimension of the 1st variable must be 2*npts*n_cam!');
Xmat = NaN(3,npts);
if n_cam<2,
    disp('WARNING: It takes at least 2 cameras to solve the 3D structure!');
    return;
end;

threshold = 1e-5;   %  threshold to assume the singularity of A
if nargin<9,
    MaxIter = 0;        % do not refine (default)
    if nargin < 8,
        alphavec = zeros(1,n_cam);
        if nargin < 7,
            kmat = zeros(5,n_cam);
        end;
    end;
end;

assert(isequal(size(Qmat),[4,n_cam]),'The 2nd variable is assumed to be quaternion matrix!');
assert(isequal(size(Tmat),[3,n_cam]),'The 3rd variable is assumed to be translation matrix!');
assert(isequal(abs(handvec),ones(1,n_cam)),'The 8th variable is assumed to behandness vector!');
[m, n] = size(fmat);
if m==2 && n==1,
    fmat = fmat(:,ones(1,n_cam));
else
    assert(m==2 && n==n_cam,'The 4th variable is assumed to be focal length matrix!');
end;

[m, n] = size(cmat);
if m==2 && n==1,
    cmat = cmat(:,ones(1,n_cam));
else
    assert(m==2 && n==n_cam,'The 5th variable is assumed to be principal point matrix!');
end;

[m, n] = size(kmat);
if m==5 && n==1,
    kmat = kmat(:,ones(1,n_cam));
else
    assert(m==5 && n==n_cam,'The 6th variable is assumed to be distortion matrix!');
end;

[m, n] = size(alphavec);
if m==1 && n==1,
    alphavec = alphavec(ones(1,n_cam));
else
    assert(m==1 && n==n_cam,'The 7th variable is assumed to be aspect ratio vector!');
end;

% logical indices of active 2d feature points
% xbody dimension (2, npts, n_cam); indmat dimension (npts, n_cam)
indmat = squeeze(all(~isnan(xbody),1));
xndata = NaN(2,npts,n_cam);
% compute projection matrix (RT):  xn=[R,T]*X
% every page is a projection matrix [R,T]
RT = zeros(3,4,n_cam);
for pp = 1:n_cam,
    % logical indices of active 2d points in each view
    indx = indmat(:,pp);
    % Compute the normalized coordinates:
    xndata(:,indx,pp) = normalize_pixel(xbody(:,indx,pp),fmat(:,pp),cmat(:,pp),kmat(:,pp),alphavec(pp));
    R = trans_quat_mat(Qmat(:,pp));
    if handvec(pp)~=1,
        R(:,3) = -R(:,3);
    end;
    RT(:,:,pp) = [R, Tmat(:,pp)];
end;

% It takes at least 2 cameras to solve the 3D structure!
% find points which can be seen in more than 2 cameras
sumind = sum(indmat,2);
ind_active_pts = find(sumind>1)';  % need row vector for loops

goodpts = false(1,npts);
for kk = ind_active_pts,
    nview = sumind(kk);
    % logical indices of active views of each feature point
    indx = indmat(kk,:);
    xn = xndata(:,kk,indx);
    rt = RT(:,:,indx);
    % The DLT method is applied here!!
    A = zeros(2,4,nview);
    A(1,:,:) = repmat(xn(1,:,:),1,4).*rt(3,:,:)-rt(1,:,:);
    A(2,:,:) = repmat(xn(2,:,:),1,4).*rt(3,:,:)-rt(2,:,:);
    A = reshape(permute(A,[2,1,3]),4,2*nview)';
    % A*X=0
    [~,s,V] = svd(A,0);
    goodpts(kk) = s(2,2)/s(1,1) >threshold && s(3,3)/s(2,2) >threshold;
    X = V(:,4);
    Xmat(:,kk) = X(1:3)/X(4);
end;

if ~MaxIter,
    return;
end;

%%%%%%%%%%%%%%%%%%%%%  LM refinement  %%%%%%%%%%%%%%%%%%%%
% if MaxIter~=0, start optimization
bigeps = eps*1e6;

Nact = length(ind_active_pts);      % active points
Nact3 = Nact*3;
Nactncam = Nact*n_cam;
% indices of good feature points
indmat = indmat(ind_active_pts,:);
indvec = reshape(indmat,1,Nactncam);
xbody = xbody(:,ind_active_pts,:);
xndata = xbody;
Xactive  = Xmat(:,ind_active_pts);
param = reshape(Xactive,Nact3,1);
for pp = 1:n_cam,
    xndata(:,:,pp) = project_points_mirror(Xactive,Qmat(:,pp),Tmat(:,pp),handvec(pp), ...
        fmat(:,pp),cmat(:,pp),kmat(:,pp),alphavec(pp));
end;
ex = reshape(xbody-xndata, 2,Nactncam);
ex = reshape(ex(:,indvec),[],1);
ex2 = dot(ex,ex);
lamda = 0.001; % set an initial value of the damping factor for the LM
updateJ = 1;

for iter = 1:MaxIter,
%     fprintf(1,'%d...',iter);
    if updateJ,
        % load parameters
        Xactive = reshape(param,3,Nact);
        % compute the approximated Hessian matrix
        JJ2 = sparse([],[],[],Nact3,Nact3,Nact3*3);     % JJ2 = JJ'*JJ;
        JJTex = zeros(Nact3,1);                                   % JJTex = JJ'*ex
        for pp = 1:n_cam,
            [xn,~,~,~,~,~,~,dxdX] = project_points_mirror(Xactive,Qmat(:,pp),Tmat(:,pp), ...
                handvec(pp),fmat(:,pp),cmat(:,pp),kmat(:,pp),alphavec(pp));
            indx = indmat(:,pp);
            dxdX = dxdX(reshape(repmat(indx',2,1),[],1), :);
            JJ2 = JJ2 + dxdX'*dxdX;
            ex = xbody(:,indx,pp)-xn(:,indx);
            JJTex = JJTex+dxdX'*ex(:);
        end;
    end;
    H_lm = JJ2 + diag(lamda*diag(JJ2));  % JJ2 + lamda*speye(Nact3);
    param_innov = H_lm\JJTex;
    param_up = param + param_innov;
    Xactive = reshape(param_up,3,Nact);
    for pp = 1:n_cam,
        xndata(:,:,pp) = project_points_mirror(Xactive,Qmat(:,pp),Tmat(:,pp),handvec(pp), ...
            fmat(:,pp),cmat(:,pp),kmat(:,pp),alphavec(pp));
    end;
    ex = reshape(xbody-xndata, 2,Nactncam);
    ex = reshape(ex(:,indvec),[],1);
    ex2_up = dot(ex,ex);
    if ex2_up < ex2,
        lamda = lamda/10;
        param = param_up;
        ex2 = ex2_up;
        updateJ=1;
    else
        lamda = lamda*10;
        updateJ=0;
    end;
    if norm(JJTex) < bigeps,                  % grad = JJTex;
        break;
    end;
end;

Xactive = reshape(param,3,Nact);
Xmat(:,ind_active_pts) = Xactive;
return;




%% test
m = 2;       % dimension of 3d points
npts = 20;
ncam = 1+randi(10);
onecam = 1;
checkoutpic = 1;
maxfov = 100;        % unit: degree

X = 100*randn(m,npts);      % unit: mm
%  Xcenter = mean(X,2); Center of X is set approximately at the origin of Xw frame
radius = 500+1000*rand(1,ncam);
% Xc=R*Xw+T, T is set approximately Z axis of Xc frame,
% so camera is looking at the origin of Xw frame.
T = [10*randn(2,ncam); radius];
% remove bad projection by set rotation angle less than pi/3
om = randn(3,ncam);
om = om./(repmat(sqrt(sum(om.*om,1)),3,1)).*repmat(rand(1,ncam),3,1)*pi/3;
q = trans_quat_axis(om);
hand = sign(randn(1,ncam));

% [x,dxdq,dxdT,dxdf,dxdc,dxdk,dxdalpha,dxdX] = project_points_mirror(X,q,T,hand,f,c,k,alpha);
xxx = NaN(2,npts,ncam);
if onecam,
    nx = 500+randi(500);
    ny = 500+randi(500);
    fov_angle = 5+randi(maxfov);
    f = ([nx;ny]/2)/tan(pi*fov_angle/360);
    c = ([nx;ny]-1)/2;
    div = 2.^(1:5)+5;
    k = randn(5,1)./div';
    alpha = 0.01*randn(1,1);
    for i = 1:ncam,
        x = project_points_mirror(X,q(:,i),T(:,i),hand(i),f,c,k,alpha);
        if checkoutpic,
            outind = x(1,:)>nx-0.5 | x(1,:)<-0.5 | x(2,:)>ny-0.5 | x(2,:)<-0.5;
            x(:,outind) = NaN;
        end;
        xxx(:,:,i) = x;
    end;
else
    imageXY = 500+randi(500,2,ncam);
    fov_angle = 5+randi(maxfov,1,ncam);
    f = (imageXY/2)./repmat(tan(pi*fov_angle/360),2,1);
    c = (imageXY-1)/2;
    div = 2.^(repmat((1:5)',1,ncam))+5;
    k = randn(5,ncam)./div;
    alpha = 0.01*randn(1,ncam);
    for i = 1:ncam,
        x = project_points_mirror(X,q(:,i),T(:,i),hand(i),f(:,i),c(:,i),k(:,i),alpha(i));
        if checkoutpic,
            nx = imageXY(1,i);
            ny = imageXY(2,i);
            outind = x(1,:)>nx-0.5 | x(1,:)<-0.5 | x(2,:)>ny-0.5 | x(2,:)<-0.5;
            x(:,outind) = NaN;
        end;
        xxx(:,:,i) = x;
    end;
end;

if m==2,
    X=[X;zeros(1,npts)];
end;
[X1,goodpts] = compute_structure(xxx,q,T,hand,f,c,k,alpha);
[X2,goodpts2] = compute_structure(xxx,q,T,hand,f,c,k,alpha,20);   % after iteration
fprintf(1,'\n');
if checkoutpic,
    indx = all(~isnan(X1),1);
    err = X1(:,indx)-X(:,indx);
    err_init = norm(err(:))
    err = X2(:,indx)-X(:,indx);
    err_lm = norm(err(:))
else
    err = X1-X;
    err_init = norm(err(:))
    err = X2-X;
    err_lm = norm(err(:))
end;

