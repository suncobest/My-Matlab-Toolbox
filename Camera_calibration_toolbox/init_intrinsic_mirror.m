%init_intrinsic_mirror: revised from init_intrinsic_param
%
%Initialization of the intrinsic parameters.
%Runs as a script.
%
%INPUT: x_cell: Feature locations on the images
%       X_cell: Corresponding grid coordinates
%
%OUTPUT: fc: Camera focal length
%        cc: Principal point coordinates
%	     kc: Distortion coefficients
%        alpha_c: skew coefficient
%        KK: The camera matrix (containing fc, cc and alpha_c)
%
%Method: Computes the planar homographies H_cell and computes
%        the focal length fc from orthogonal vanishing points constraint.
%        The principal point cc is assumed at the center of the image.
%        Assumes no image distortion (kc = [0;0;0;0])
%
%Note: The row vector active_images consists of zeros and ones. To deactivate an image, set the
%      corresponding entry in the active_images vector to zero.
%
%
%Important function called within that program:
%
%compute_homography.m: Computes the planar homography between points on the grid in 3D, and the image plane.
%
%
%VERY IMPORTANT: This function works only with 2D rigs and projective
%cameras. This code is not in working for affine cameras.
%对于射影相机有：x=PX，其中P=sKR[I|-C]=sK[r1,r2,r3,t]。对于仿射相机，P=sK[R(1:2,:);0 0 0],t]
%棋盘格与像平面的单应矩阵变成仿射单应，由射影变换H=sK[r1,r2,t]变为仿射变换H=sK[R(1:2,1:2);0 0],t]
%射影相机成为仿射相机的过程：长焦距或者窄视野以及远距离拍摄

if ~exist('x_cell','var'),
    fprintf(1,['\nThere is no data required for intrinsic calibration!']);
    click_track_mirror;
end;

if ~exist('nx','var') || ~exist('ny','var'),
    fprintf(1,'WARNING: No image size (nx,ny) available. Please set the image size manually!\n');
    nx = input('Width: nx ([] = 640) = ');
    if isempty(nx),
        nx = 640;
    end;
    ny = input('Height: ny ([] = 480) = ');
    if isempty(ny),
        ny = 480;
    end;
end;

if ~exist('est_aspect_ratio','var'),
    est_aspect_ratio = 1;
end;

if ~exist('est_fc','var'),
    est_fc = [1;1];
end;

% Initialize the homographies:
n_view = n_ima * n_cam;
H_cell = cell(1,n_view);
for kk = 1:n_view,
    x_kk = x_cell{kk};
    X_kk = X_cell{kk};
    if isempty(x_kk),
        if active_imgviews(kk),
            fprintf(1,'Warning: Data is missing! Set active_imgviews(%d)=0;\n',kk);
        end;
        active_imgviews(kk) = 0;
    end;
    if active_imgviews(kk),
        H_cell{kk}= compute_homography_lm(x_kk,X_kk(1:2,:));
    end;
end;

kk = find(active_images ~= any(active_imgviews,1));
if ~isempty(kk),
    fprintf(1,'WARNING: Cannot calibrate all views of image %d.\n',kk);
    fprintf(1,'Set active_images(%d)=0;\n',kk);
end;
active_images = any(active_imgviews);
ind_active_views = find(active_imgviews(:))';

fprintf(1,'\nInitialization of the intrinsic parameters - Number of image views: %d\n',length(ind_active_views));
% initial guess of distortion:
k_init = [0;0;0;0;0]; % initialize to zero (no distortion)
% Compute K form scene and internal constraints
% n1'*IAC*h2=0; h1'*IAC*h1=h2'*IAC*h2; IAC=inv(K')*inv(K) 
% IAC=[w1,w2,w4;w2,w3,w5;w4,w5,w6], IACv = [w1,w2,w3,w4,w5,w6]'

% build the matrix:
A = [];
for kk = ind_active_views,
    Hkk = H_cell{kk};     % H = [h1,h2,h3]=sK(r1,r2,t)
    h1 = Hkk(:,1);
    h2 = Hkk(:,2);
    A12 = [h1(1)*h2(1), h1(1)*h2(2)+h1(2)*h2(1), h1(2)*h2(2), h1(1)*h2(3)+h1(3)*h2(1), h1(2)*h2(3)+h1(3)*h2(2), h1(3)*h2(3)];
    A11 = [h1(1)*h1(1), h1(1)*h1(2)*2, h1(2)*h1(2), h1(1)*h1(3)*2, h1(2)*h1(3)*2, h1(3)*h1(3)];
    A22 = [h2(1)*h2(1), h2(1)*h2(2)*2, h2(2)*h2(2), h2(1)*h2(3)*2, h2(2)*h2(3)*2, h2(3)*h2(3)];
    A = [A; A12; A11-A22];
end;

% A*IACv=0
A = A'*A;
[~,~,V] = svd(A);
IACv = V(:,6);  % 最小奇异值对应的右奇异矢量，且norm(hh)=1
if IACv(6)<0,
    IACv = -IACv;
end;
IAC = [IACv(1),IACv(2),IACv(4); IACv(2),IACv(3),IACv(5); IACv(4),IACv(5),IACv(6)];   % IAC为绝对二次曲线的像

% Cholesky factorization
inv_KK = chol(IAC);
KK = inv(inv_KK);
KK = KK/KK(end);
inv_KK = inv(KK);

% clear A Hkk h1 h2 A12 A11 A22 U S V;

cc = KK(7:8)';
fc = [KK(1); KK(5)];
if ~est_aspect_ratio && all(est_fc),
    fc(1) = (fc(1)+fc(2))/2;
    fc(2) = fc(1);
end;
kc = k_init;
alpha_c = KK(4)/KK(1);

if ((cc(1)<-.5) || (cc(1)>nx-.5) || (cc(2)<-.5) || (cc(2)>ny-.5)),
    % It takes a large shift for View camera lens to have this result, and it''s almost impossible for an ordinary camera!
    fprintf(1,['Warning: cc = [ %3.5f   %3.5f ].\nIt appears that the principal point is out of the picture!\n' ...
        'Your lens must shift out of the picture!\n'],cc);
    flag = input('Set principal point to the image center or not? ([]=no, other=yes) ','s');
    if ~isempty(flag),
        fprintf(1,'Setting principal point to the center of the image!\n');
        cc = [nx;ny]/2 - 0.5; % initialize at the center of the image
    end;
end;

% KK = [f_init(1) alpha_init*f_init(1) c_init(1);0 f_init(2) c_init(2); 0 0 1];
% alpha_c=-ctg(theta), theta is the angle between the pixel axes x and y （左上角）
% angle of pixel = 180 - theta（左下角）

fprintf(1,'\n\nCalibration parameters after initialization:\n\n');
fprintf(1,'Focal Length:          fc = [%3.5f, %3.5f]\n',fc);
fprintf(1,'Principal point:       cc = [%3.5f, %3.5f]\n',cc);
fprintf(1,'Skew:             alpha_c =  %3.5f  => Skew angle = %3.5f degrees\n',alpha_c,90-atan(alpha_c)*180/pi);
fprintf(1,'Distortion:              kc = [%3.5f, %3.5f, %3.5f, %3.5f, %3.5f]\n',kc);
