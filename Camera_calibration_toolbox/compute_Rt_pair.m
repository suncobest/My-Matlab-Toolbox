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

assert(isequal(size(omcmat),[3,n_cam]),'The 2nd variable is assumed to be axis angle matrix!');
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

[x1, T1] = normalise2dpts(x1);
[x2, T2] = normalise2dpts(x2);

% Build the constraint matrix
A = [x2(1,:)'.*x1(1,:)',  x2(1,:)'.*x1(2,:)',  x2(1,:)', ...
     x2(2,:)'.*x1(1,:)',  x2(2,:)'.*x1(2,:)',  x2(2,:)', ...
     x1(1,:)',             x1(2,:)',           ones(npts,1)];

[U,D,V] = svd(A,0);

% Extract fundamental matrix from the column of V corresponding to
% smallest singular value.
F = reshape(V(:,9),3,3)';

% Enforce constraint that fundamental matrix has rank 2 by performing
% a svd and then reconstructing with the two largest singular values.
[U,D,V] = svd(F,0);
F = U*diag([D(1,1) D(2,2) 0])*V';

% Denormalise
F = T2'*F*T1;

if nargout > 1,  	% Solve for epipoles
    [U,D,V] = svd(F,0);
    if V(3,3)>eps,
        e1 = V(:,3)/V(3,3);
    end;
    if U(3,3)>eps,
        e2 = U(:,3)/U(3,3);
    end;
    if nargout==2,
        e1 = [e1,e2];
    end;
end;

return;



%% Test
