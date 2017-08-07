function [Omc,Tc,Emat] = compute_Rt_pair(xpair,fmat,cmat,kmat,alphavec,hand,MaxIter)
% COMPUTE_RT_PAIR - computes essential matrix from corresponding points of two cameras
% and decompose it into R and t. Where
%         Xc2 = R*hand*Xc1+t.
% Function computes the essential matrix from 8 or more matching points in
% a stereo pair of images.  The normalised 8 point algorithm given by
% Hartley and Zisserman p265 is used.
%
% Arguments:
%          xpair: corresponding image points of two cameras: 2*npts*2 or 4*npts
%          fmat: Camera focal length of every camera: 2*2
%          cmat: Principal point coordinates of every Camera: 2*2
%          kmat: Distortion coefficients of every camera: 5*2
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

if nargin < 7,
    MaxIter = 0;        % do not refine (default)
    if nargin < 6,
        hand = 1;
        if nargin < 5,
            alphavec = zeros(1,2);
            if nargin < 4,
                kmat = zeros(5,2);
            end;
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

x1n = [normalize_pixel(x1,fmat(:,1),cmat(:,1),kc(:,1),alphavec(1)); ones(1,npts)];
x2n = [normalize_pixel(x2,fmat(:,2),cmat(:,2),kc(:,2),alphavec(2)); ones(1,npts)];

% Build the constraint matrix
A = [x2n(1,:)'.*x1n(1,:)',  x2n(1,:)'.*x1n(2,:)',  x2n(1,:)', ...
     x2n(2,:)'.*x1n(1,:)',  x2n(2,:)'.*x1n(2,:)',  x2n(2,:)', ...
     x1n(1,:)',             x1n(2,:)',           ones(npts,1)];

[U,D,V] = svd(A,0);

% Extract fundamental matrix from the column of V corresponding to
% smallest singular value.
E = reshape(V(:,9),3,3)';

% Enforce constraint that fundamental matrix has rank 2 by performing
% a svd and then reconstructing with the two largest singular values.
[U,D,V] = svd(E,0);
E = U*diag([1, 1, 0])*V';


return;



%% Test
