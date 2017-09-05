%init_intrinsic_mirror2: revised from init_intrinsic_param
%
%Initialization of the intrinsic parameters.
%Runs as a script.
%
%INPUT: x_1,x_2,x_3,...: Feature locations on the images
%       X_1,X_2,X_3,...: Corresponding grid coordinates
%
%OUTPUT: fc: Camera focal length
%        cc: Principal point coordinates
%	     kc: Distortion coefficients
%        alpha_c: skew coefficient
%        KK: The camera matrix (containing fc, cc and alpha_c)
%
%Method: Computes the planar homographies H_1, H_2, H_3, ... and computes
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
%VERY IMPORTANT: This function works only with 2D rigs.
%In the future, a more general function will be there (working with 3D rigs as well).


if ~exist('est_aspect_ratio','var'),
    est_aspect_ratio = 1;
end;

check_active_images;
ind_active_views = find(active_imgviews(:))';

if ~exist(['x_' num2str(ind_active_views(1))],'var'),
    click_calib_mirror;
end;


fprintf(1,'\nInitialization of the intrinsic parameters - Number of image views: %d\n',length(ind_active_views));


% Initialize the homographies:

for kk = 1:n_ima,
    for pp = 1:n_cam,
        eval(['x_kk = x_' num2str((kk-1)*n_cam+pp) ';']);
        eval(['X_kk = X_' num2str((kk-1)*n_cam+pp) ';']);
        if (isnan(x_kk(1,1))),
            if active_imgviews(pp,kk),
                fprintf(1,'Warning: Cannot calibrate with view %d of image %d.\n',pp,kk);
                fprintf(1,'         Set active_imgviews(%d,%d)=0;\n',pp,kk);
                active_imgviews(pp,kk) = 0;
            end;
        end;
        if active_imgviews(pp,kk),
            eval(['H_' num2str((kk-1)*n_cam+pp) ' = compute_homography_lm(x_kk,X_kk(1:2,:));']);
        else
            eval(['H_' num2str((kk-1)*n_cam+pp) ' = NaN*ones(3,3);']);
        end;
    end;

    if active_images(kk)~= any(active_imgviews(:,kk)),
        fprintf(1,'WARNING: Cannot calibrate all views of image %d.\n',kk);
        fprintf(1,'         Set active_images(%d)=0;\n',kk);
    end;
end;

active_images = any(active_imgviews);
check_active_images;

% initial guess for principal point and distortion:

if ~exist('nx','var') || ~exist('ny','var'), [ny,nx] = size(I); end;

k_init = [0;0;0;0;0]; % initialize to zero (no distortion)
c_init = [nx;ny]/2 - 0.5; % initialize at the center of the image

A = [];
b = [];

% matrix that subtract the principal point:
Sub_cc = [1 0 -c_init(1);0 1 -c_init(2);0 0 1];
n_view = n_ima * n_cam;

for kk=1:n_view,

    if active_imgviews(kk),

        eval(['Hkk = H_' num2str(kk) ';']);

        Hkk = Sub_cc * Hkk;
        %将图像坐标原点平移到主点，不考虑倾斜alpha_c，则有Hkk = diag(fc1,fc2,1)*[r1,r2,t]
        %此时内参数K=diag(fc1,fc2,1)；绝对二次曲线的像（IAC）为omega=inv(K')*inv(K) = diag(1/fc1^2,1/fc2^2,1)
        % 根据r1'*r2=0和r1'*r1=r2'*r2，所以有[a1,b1,c1]*diag(1/fc1^2,1/fc2^2,1)*[a2;b2;c2]=0；
        % [a1,b1,c1]*diag(1/fc1^2,1/fc2^2,1)*[a1;b1;c1] = [a2,b2,c2]*diag(1/fc1^2,1/fc2^2,1)*[a2;b2;c2]
        % 所以[a1*a2, b1*b2, c1*c2]*[1/fc1^2;1/fc2^2;1]=0, [a1*a1-a2*a2, b1*b1-b2*b2, c1*c1-c2*c2]*[1/fc1^2;1/fc2^2;1]=0
        % A_kk*[1/fc1^2;1/fc2^2]=b_kk，增广矩阵A和b也满足A*[1/fc1^2;1/fc2^2]=b

        % [a1*a1-a2*a2, b1*b1-b2*b2, c1*c1-c2*c2]*[1/fc1^2;1/fc2^2;1]=[a1+a2,b1+b2,c1+c2]*diag(1/fc1^2,1/fc2^2,1)*[a1-a2;b1-b2;c1-c2]=0；
        % 则h1+h2和h1-h2为共轭消失点
        % 这里的算法等价于init_intrinsic_param.m

        V_hori_pix = Hkk(:,1);    % h1
        V_vert_pix = Hkk(:,2);    % h2

        a1 = V_hori_pix(1);
        b1 = V_hori_pix(2);
        c1 = V_hori_pix(3);

        a2 = V_vert_pix(1);
        b2 = V_vert_pix(2);
        c2 = V_vert_pix(3);

        A_kk = [a1*a2  b1*b2;
            a1^2-a2^2  b1^2-b2^2];

        b_kk = -[c1*c2;c1^2-c2^2];


        A = [A;A_kk];
        b = [b;b_kk];

    end;
end;

% keyboard;

if est_aspect_ratio,
    % Use a two focals estimate:
    f_init = sqrt(abs(1./((A'*A)\(A'*b)))); % if using a two-focal model for initial guess，其中f_init = [fc1;fc2]
else
    % Use a single focal estimate:
    A1 = sum(A,2);
    f_init = sqrt((A1'*A1)/(A1'*b))* ones(2,1); % if single focal length model is used
end;

alpha_init = 0;

% Global calibration matrix (initial guess):

KK = [f_init(1) alpha_init*f_init(1) c_init(1);0 f_init(2) c_init(2); 0 0 1];
inv_KK = inv(KK);


cc = c_init;
fc = f_init;
kc = k_init;
alpha_c = alpha_init;

% alpha_c=-ctg(theta), theta is the angle between the pixel axes x and y （左上角）
% angle of pixel = 180 - theta（左下角）

fprintf(1,'\n\nCalibration parameters after initialization:\n\n');
fprintf(1,'Focal Length:          fc = [ %3.5f   %3.5f ]\n',fc);
fprintf(1,'Principal point:       cc = [ %3.5f   %3.5f ]\n',cc);
fprintf(1,'Skew:             alpha_c = [ %3.5f ]   => angle of pixel = %3.5f degrees\n',alpha_c,90 - atan(alpha_c)*180/pi);
fprintf(1,'Distortion:            kc = [ %3.5f   %3.5f   %3.5f   %3.5f   %5.5f ]\n',kc);
